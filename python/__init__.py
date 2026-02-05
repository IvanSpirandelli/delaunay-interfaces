"""
DelaunayInterfaces Python Bindings

This module provides Python bindings for the DelaunayInterfaces C++ library,
which computes interface surfaces from multicolored point clouds using
Delaunay/alpha complexes and barycentric subdivision.
"""

from .delaunay_interfaces import (
    InterfaceGenerator,
    InterfaceSurface,
    get_barycentric_subdivision_and_filtration,
    __version__
)

__all__ = [
    'InterfaceGenerator',
    'InterfaceSurface',
    'get_barycentric_subdivision_and_filtration',
]
