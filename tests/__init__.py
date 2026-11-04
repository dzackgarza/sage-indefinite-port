"""The port test tree: a regular package so the acceptance suite imports its own modules.

An installed third-party package in the pinned image owns a regular ``tests`` package;
a namespace-style directory of the same name loses to it on ``sys.path``.
"""
