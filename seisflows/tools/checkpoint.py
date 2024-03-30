"""
A simple checkpointing system that is used for failure tolerance and seamless
restarts during workflows. Keeps track of both function names and system task
IDs so that functions which partially finish do not rerun already completed
tasks when resuming
"""
import numpy as np


class Checkpoint:
    """
    Simple checkpointing system for completed tasks and sub-tasks
    """
