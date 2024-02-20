"""
Preprocessing tools for the ambient noise adjoint tomography Noise
Inversion workflow module. Tools here are used for waveform rotation
required for horizontal component sensitivity kernels (RR and TT).
"""
import numpy as np


def rotate_ne_trace_to_rt(tr_ee, tr_ne, tr_en, tr_nn, theta, theta_p):
    """
    Used during ambient noise adjoint tomography (workflow: noise inversion) to
    rotate N and E component Synthetic Greens Functions (SGF) to R and T
    components so that they can be directly compared to RR and TT
    Empirical Greens Functions (EGF)

    Naming convention: AB (A=force source direction, B=waveform component)
    Order of input follows Wang et al. 2019 Eqs. 9 and 10

    :type tr_ee: obspy.core.trace.Trace
    :param tr_ee: E component FORCE, E component SGF
    :type tr_ne: obspy.core.trace.Trace
    :param tr_ne: N component FORCE, E component SGF
    :type tr_en: obspy.core.trace.Trace
    :param tr_en: E component FORCE, N component SGF
    :type tr_nn: obspy.core.trace.Trace
    :param tr_nn: N component FORCE, N component SGF
    :type theta: float
    :param theta: azimuth between source station and receiver station in units
        of radians (see Wang et al. 2018 Fig. 1)
    :type theta_p: float
    :param theta_p: theta prime, 180 degrees from the backazimuth of the source
        station and receiver station in units radians. theta != theta_prime
        for a spherical Earth, but they will be close.
    """
    tr_tt = tr_nn.copy()
    tr_rr = tr_nn.copy()

    # TT rotation from Wang et al. (2019) Eq. 9
    tr_tt.data = (+ 1 * np.cos(theta) * np.cos(theta_p) * tr_ee.data
                  - 1 * np.cos(theta) * np.sin(theta_p) * tr_ne.data
                  - 1 * np.sin(theta) * np.cos(theta_p) * tr_en.data
                  + 1 * np.sin(theta) * np.sin(theta_p) * tr_nn.data
                  )
    tr_tt.stats.component = "T"

    # RR rotation from Wang et al. (2019) Eq. 10
    tr_rr.data = (+ 1 * np.sin(theta) * np.sin(theta_p) * tr_ee.data
                  + 1 * np.sin(theta) * np.cos(theta_p) * tr_ne.data
                  + 1 * np.cos(theta) * np.sin(theta_p) * tr_en.data
                  + 1 * np.cos(theta) * np.cos(theta_p) * tr_nn.data
                  )
    tr_rr.stats.component = "R"

    return tr_rr, tr_tt


def rotate_rt_adjsrc_to_ne(tr, theta, theta_p):
    """
    Rotate RR and TT adjoint sources to N and E components for use in
    N and E component adjoint simulations in ambient noise inversion workflow.

    Paraphrasing Wang et al. (2019) 2.2.2

    Kernels are generated from the time convolution of a forward field due to a
    point force in the [east or north] direction, and an adjoint field obtained
    by injecting [f_TT or f_RR (adjoint source)] in the east and north
    directions at the receiver, with their amplitude modulated by the rotation
    matrix.

    :type tr: obspy.core.trace.Trace
    :param tr: input trace that is either R or T component. Component stat
        of the Trace is required (Trace.stats.component)
    :type theta: float
    :param theta: azimuth between source station and receiver station in units
        of radians
    :type theta_p: float
    :param theta_p: theta prime, 180 degrees from the backazimuth of the source
        station and receiver station in units radians. theta != theta_prime
        for a spherical Earth, but they will be close.
    :rtype: tuple of obspy.core.trace.Trace
    :return: (EE, NE, EN, NN) component adjoint sources that have been rotated
        based on the component of `tr_in`
    """
    component = tr.stats.component
    assert(component in ["R", "T"]), \
        "input adjoint source must have component 'R' or 'T'"

    # Naming convention: AB (A=force source direction, B=injection direction)
    # Copy the input trace to modify for the output traces
    tr_ee = tr.copy()
    tr_en = tr.copy()
    tr_ne = tr.copy()
    tr_nn = tr.copy()

    # TT rotation from Wang et al. (2019) Eq. 16
    # Note that, confusingly, the rotation matrix notation order changes for
    # w.r.t synthetics rotation ('en' switches with 'ne'), because the
    # summation order changes for adjoint source rotation
    if component == "T":
        tr_ee.data = + 1 * np.cos(theta) * np.cos(theta_p) * tr.data  # K_ET
        tr_en.data = - 1 * np.cos(theta) * np.sin(theta_p) * tr.data  # K_ET
        tr_ne.data = - 1 * np.sin(theta) * np.cos(theta_p) * tr.data  # K_NT
        tr_nn.data = + 1 * np.sin(theta) * np.sin(theta_p) * tr.data  # K_NT
    # TT rotation from Wang et al. (2019) Eq. 18
    elif component == "R":
        tr_ee.data = + 1 * np.sin(theta) * np.sin(theta_p) * tr.data  # K_ER
        tr_en.data = + 1 * np.sin(theta) * np.cos(theta_p) * tr.data  # K_ER
        tr_ne.data = + 1 * np.cos(theta) * np.sin(theta_p) * tr.data  # K_NR
        tr_nn.data = + 1 * np.cos(theta) * np.cos(theta_p) * tr.data  # K_NR

    # Return in order that all other functions are using, despite the rotation
    # matrix having a different order
    return tr_ee, tr_ne, tr_en, tr_nn
