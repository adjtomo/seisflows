#!/usr/bin/env python3
"""
Misfit functions used by the 'default' preprocess class use to quantify 
differences between data and synthetics. 

All functions defined have four required positional arguments

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
"""
import numpy as np
from scipy.signal import hilbert as analytic


def waveform(syn, obs, nt, dt, *args, **kwargs):
    """
    Direct waveform differencing

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    wrsd = syn - obs

    return np.sqrt(np.sum(wrsd * wrsd * dt))


def envelope(syn, obs, nt, dt, *args, **kwargs):
    """
    Waveform envelope difference from Yuan et al. 2015 Eq. 9

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    env_syn = abs(analytic(syn))
    env_obs = abs(analytic(obs))

    # Residual of envelopes
    env_rsd = env_syn - env_obs

    return np.sqrt(np.sum(env_rsd * env_rsd * dt))


def instantaneous_phase(syn, obs, nt, dt, *args, **kwargs):
    """
    Instantaneous phase difference from Bozdag et al. 2011

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    r = np.real(analytic(syn))
    i = np.imag(analytic(syn))
    phi_syn = np.arctan2(i, r)

    r = np.real(analytic(obs))
    i = np.imag(analytic(obs))
    phi_obs = np.arctan2(i, r)

    phi_rsd = phi_syn - phi_obs

    return np.sqrt(np.sum(phi_rsd * phi_rsd * dt))


def traveltime(syn, obs, nt, dt, *args, **kwargs):
    """
    Cross-correlation traveltime 

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    cc = abs(np.convolve(obs, np.flipud(syn)))

    return (np.argmax(cc) - nt + 1) * dt


def traveltime_inexact(syn, obs, nt, dt, *args, **kwargs):
    """
    A faster cc traveltime function but possibly innacurate

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    it = np.argmax(syn)
    jt = np.argmax(obs)

    return (jt - it) * dt


def amplitude(syn, obs, nt, dt, *args, **kwargs):
    """
    Cross-correlation amplitude difference

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    ioff = (np.argmax(cc) - nt + 1) * dt

    if ioff <= 0:
        wrsd = syn[ioff:] - obs[:-ioff]
    else:
        wrsd = syn[:-ioff] - obs[ioff:]

    return np.sqrt(np.sum(wrsd * wrsd * dt))


def envelope2(syn, obs, nt, dt, *args, **kwargs):
    """
    Envelope amplitude ratio from Yuan et al. 2015 Eq. B-1

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    env_syn = abs(analytic(syn))
    env_obs = abs(analytic(obs))

    raise NotImplementedError


def envelope3(syn, obs, nt, dt, eps=0., *args, **kwargs):
    """
    Envelope cross-correlation lag from Yuan et al. 2015, Eq. B-4

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    env_syn = abs(analytic(syn))
    env_obs = abs(analytic(obs))

    return Traveltime(env_syn, env_obs, nt, dt)


def instantaneous_phase2(syn, obs, nt, dt, eps=0., *args, **kwargs):
    """
    Alterative instantaneous phase function

    :type syn: np.array
    :param syn: synthetic data array
    :type obs: np.array
    :param obs: observed data array
    :type nt: int
    :param nt: number of time steps in the data array
    :type dt: float
    :param dt: time step in sec
    """
    env_syn = abs(analytic(syn))
    env_obs = abs(analytic(obs))

    env_syn1 = env_syn + eps * max(env_syn)
    env_obs1 = env_obs + eps * max(env_obs)

    diff = (syn / env_syn1) - (obs / env_obs1)

    return np.sqrt(np.sum(diff * diff * dt))


def displacement(*args, **kwargs):
    return Exception("This function can only used for migration.")


def velocity(*args, **kwargs):
    return Exception("This function can only used for migration.")


def acceleration(*args, **kwargs):
    return Exception("This function can only used for migration.")

