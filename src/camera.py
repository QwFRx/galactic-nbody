"""Orbiting perspective camera shared by the live view and the renderer."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Camera:
    """Simple yaw-pitch-distance camera looking at a target point."""

    yaw: float = 35.0
    pitch: float = 25.0
    distance: float = 400.0
    target: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=np.float64))
    fov: float = 50.0
    min_distance_scale: float = 0.15
    max_distance_scale: float = 8.0
    _extent: float = 300.0

    def reset(self, target: np.ndarray, extent: float) -> None:
        """Place the camera at a pleasant three-quarter view of the scene."""
        self.target = np.asarray(target, dtype=np.float64).copy()
        self._extent = float(extent)
        self.distance = 1.35 * float(extent)
        self.yaw = 35.0
        self.pitch = 25.0

    def eye(self) -> np.ndarray:
        cy, sy = np.cos(np.radians(self.yaw)), np.sin(np.radians(self.yaw))
        cp, sp = np.cos(np.radians(self.pitch)), np.sin(np.radians(self.pitch))
        direction = np.array([cp * cy, cp * sy, sp], dtype=np.float64)
        return self.target + self.distance * direction

    def basis(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (forward, right, up) unit vectors."""
        forward = self.target - self.eye()
        norm = np.linalg.norm(forward)
        forward = forward / max(norm, 1.0e-12)
        world_up = np.array([0.0, 0.0, 1.0])
        right = np.cross(forward, world_up)
        if np.linalg.norm(right) < 1.0e-6:
            right = np.array([1.0, 0.0, 0.0])
        right = right / np.linalg.norm(right)
        up = np.cross(right, forward)
        return forward, right, up

    def projection(self, width: int, height: int
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        forward, right, up = self.basis()
        focal = 0.5 * height / np.tan(0.5 * np.radians(self.fov))
        return self.eye(), forward, right, focal

    def orbit(self, dyaw: float, dpitch: float) -> None:
        self.yaw += float(dyaw)
        self.pitch = float(np.clip(self.pitch + dpitch, -89.0, 89.0))

    def zoom(self, factor: float) -> None:
        lo = self.min_distance_scale * self._extent
        hi = self.max_distance_scale * self._extent
        self.distance = float(np.clip(self.distance * factor, lo, hi))

    def project(self, points: np.ndarray, width: int, height: int
                ) -> tuple[np.ndarray, np.ndarray]:
        """Project points to pixel coordinates; returns (xy, depth)."""
        eye, forward, right, focal = self.projection(width, height)
        rel = np.asarray(points, dtype=np.float64) - eye
        depth = rel @ forward
        safe = np.maximum(depth, 1.0e-9)
        sx = width * 0.5 + (rel @ right) / safe * focal
        sy = height * 0.5 - (rel @ np.cross(right, forward)) / safe * focal
        return np.stack((sx, sy), axis=1), depth
