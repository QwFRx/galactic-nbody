"""Pygame front end: orbiting camera, live panel, help and frame recording.

The simulation runs in :class:`src.sim_thread.SimulationThread`.  This module
only imports pygame when the interactive mode is actually selected, so the
headless CLI never creates a window or initialises an audio/video subsystem.
"""

from __future__ import annotations

import os
import time
from typing import Any

import numpy as np
import pygame

import config
from .camera import Camera
from .render import Renderer, save_png
from .sim_thread import SimulationThread

COLOR_MODES = ("galaxy", "component", "density", "speed")


class App:
    """Interactive viewer for either simulation scene."""

    def __init__(self, cfg: config.SimConfig,
                 thread: SimulationThread | None = None) -> None:
        self.cfg = cfg
        self.thread = thread or SimulationThread(cfg)
        self.width = cfg.win_w
        self.height = cfg.win_h
        pygame.init()
        self.screen = pygame.display.set_mode(
            (self.width, self.height), pygame.RESIZABLE)
        pygame.display.set_caption("Galactic N-Body — collision / cosmic web")
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont("dejavusans", 15)
        self.small = pygame.font.SysFont("dejavusans", 13)
        self.mono = pygame.font.SysFont("dejavusansmono", 14)
        self.big = pygame.font.SysFont("dejavusans", 22, bold=True)

        self.renderer = Renderer(cfg.render_w, cfg.render_h)
        self.hires = Renderer(1920, 1080)
        self.camera = Camera()
        self.camera_ready = False
        self.color_mode = 0
        self.show_dm = self.cfg.scene == "collision"
        self.exposure = 5.0
        self.slab = 0.35
        self.show_help = False
        self.running = True
        self.dragging = False
        self.recording = False
        self.record_index = 0
        self.record_tick = 0
        self.record_dir = os.path.join("images", "recording")
        self.last_frame: np.ndarray | None = None
        self.last_snapshot_id: int | None = None
        self.last_signature: tuple[Any, ...] | None = None
        self.status = "Initialising simulation..."
        self.status_until = time.time() + 4.0
        self.thread.start()

    # -- lifecycle ---------------------------------------------------------
    def run(self) -> None:
        """Enter the pygame event/render loop."""
        while self.running:
            self._events()
            self._update_camera()
            self._draw()
            self.clock.tick(60)
        self.thread.stop()
        pygame.quit()

    # -- input -------------------------------------------------------------
    def _events(self) -> None:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                self.running = False
            elif event.type == pygame.VIDEORESIZE:
                self.width = max(320, event.w)
                self.height = max(200, event.h)
                self.screen = pygame.display.set_mode(
                    (self.width, self.height), pygame.RESIZABLE)
            elif event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                self.dragging = True
            elif event.type == pygame.MOUSEBUTTONUP and event.button == 1:
                self.dragging = False
            elif event.type == pygame.MOUSEMOTION and self.dragging:
                self.camera.orbit(event.rel[0] * 0.32,
                                  -event.rel[1] * 0.32)
            elif event.type == pygame.MOUSEWHEEL:
                self.camera.zoom(0.90 if event.y > 0 else 1.11)
            elif event.type == pygame.KEYDOWN:
                self._key(event)

    def _key(self, event: pygame.event.Event) -> None:
        key = event.key
        if key == pygame.K_ESCAPE:
            self.running = False
        elif key == pygame.K_SPACE:
            paused = self.thread.toggle_pause()
            self._toast("Paused" if paused else "Running")
        elif key == pygame.K_h:
            self.show_help = not self.show_help
        elif key == pygame.K_c:
            self.color_mode = (self.color_mode + 1) % len(COLOR_MODES)
            self._toast(f"Colour mode: {COLOR_MODES[self.color_mode]}")
        elif key == pygame.K_d and self.cfg.scene == "collision":
            self.show_dm = not self.show_dm
            self._toast(f"Dark matter: {'on' if self.show_dm else 'off'}")
        elif key == pygame.K_r:
            self._reset_camera()
            self._toast("Camera reset")
        elif key in (pygame.K_PLUS, pygame.K_EQUALS, pygame.K_KP_PLUS):
            self.thread.set_speed(self.thread.speed * 1.5)
            self._toast(f"Speed: {self.thread.speed:.2f}x")
        elif key in (pygame.K_MINUS, pygame.K_KP_MINUS):
            self.thread.set_speed(self.thread.speed / 1.5)
            self._toast(f"Speed: {self.thread.speed:.2f}x")
        elif key == pygame.K_LEFT:
            self.camera.orbit(-4.0, 0.0)
        elif key == pygame.K_RIGHT:
            self.camera.orbit(4.0, 0.0)
        elif key == pygame.K_UP:
            self.camera.orbit(0.0, 3.0)
        elif key == pygame.K_DOWN:
            self.camera.orbit(0.0, -3.0)
        elif key == pygame.K_LEFTBRACKET and self.cfg.scene == "cosmic":
            self.slab = max(0.05, self.slab - 0.05)
            self._toast(f"Slab: {self.slab:.2f} of the box")
        elif key == pygame.K_RIGHTBRACKET and self.cfg.scene == "cosmic":
            self.slab = min(1.00, self.slab + 0.05)
            self._toast(f"Slab: {self.slab:.2f} of the box")
        elif key == pygame.K_s:
            self._save_frame()
        elif key == pygame.K_f:
            self._save_hires()
        elif key == pygame.K_v:
            self.recording = not self.recording
            if self.recording:
                os.makedirs(self.record_dir, exist_ok=True)
            self._toast("Recording: " + ("on" if self.recording else "off"))
        elif key == pygame.K_t:
            self._restart()

    def _toast(self, text: str) -> None:
        self.status = text
        self.status_until = time.time() + 3.0

    # -- simulation helpers ------------------------------------------------
    def _snapshot(self) -> dict[str, Any] | None:
        return self.thread.latest()

    def _reset_camera(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        extent = self._extent(snap)
        self.camera.reset(np.asarray(snap["center"], dtype=np.float64), extent)
        self.camera_ready = True

    def _extent(self, snap: dict[str, Any]) -> float:
        box = float(snap["box"])
        if snap["scene"] == "cosmic":
            return 0.82 * box
        return 0.28 * box

    def _update_camera(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        if not self.camera_ready:
            self._reset_camera()
        # The collision box grows adaptively; keep the camera framing while
        # the user is not actively zooming by updating only the extent.
        if snap["scene"] == "collision":
            target = np.asarray(snap["center"], dtype=np.float64)
            if np.linalg.norm(target - self.camera.target) > 1.0e-3:
                self.camera.target = target
        else:
            self.camera.target = np.asarray(snap["center"], dtype=np.float64)

    def _restart(self) -> None:
        self.thread.stop()
        self.cfg.seed += 1
        self.thread = SimulationThread(self.cfg)
        self.thread.start()
        self.camera_ready = False
        self.last_frame = None
        self.last_snapshot_id = None
        self.last_signature = None
        self._toast(f"Restarted with seed {self.cfg.seed}")

    # -- rendering ---------------------------------------------------------
    def _make_frame(self, snap: dict[str, Any]) -> np.ndarray:
        signature = (self.thread.version, self.color_mode, self.show_dm,
                     round(self.exposure, 3), round(self.slab, 3),
                     tuple(np.round(self.camera.target, 3)),
                     round(self.camera.yaw, 2), round(self.camera.pitch, 2),
                     round(self.camera.distance, 3))
        if self.last_frame is not None and signature == self.last_signature:
            return self.last_frame
        frame = self.renderer.render(
            snap, self.camera, self.color_mode, self.show_dm,
            exposure=self.exposure, slab=self.slab)
        self.last_frame = frame
        self.last_signature = signature
        if self.recording:
            self.record_tick += 1
            if self.record_tick % 2 == 0:
                os.makedirs(self.record_dir, exist_ok=True)
                path = os.path.join(self.record_dir,
                                    f"frame_{self.record_index:05d}.png")
                save_png(frame, path)
                self.record_index += 1
        return frame

    def _save_frame(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        frame = self.renderer.render(snap, self.camera, self.color_mode,
                                     self.show_dm, self.exposure, self.slab)
        path = os.path.join("images",
                            time.strftime("render_%Y%m%d_%H%M%S.png"))
        save_png(frame, path)
        self._toast(f"Saved {path}")

    def _save_hires(self) -> None:
        snap = self._snapshot()
        if snap is None:
            return
        frame = self.hires.render(snap, self.camera, self.color_mode,
                                  self.show_dm, self.exposure, self.slab)
        path = os.path.join("images",
                            time.strftime("frame_1920_%Y%m%d_%H%M%S.png"))
        save_png(frame, path)
        self._toast(f"Saved {path} (1920x1080)")

    # -- drawing -----------------------------------------------------------
    def _draw(self) -> None:
        self.screen.fill((4, 6, 12))
        snap = self._snapshot()
        if snap is not None:
            if self.status == "Initialising simulation...":
                self.status_until = 0.0
            frame = self._make_frame(snap)
            source = pygame.surfarray.make_surface(
                np.transpose(frame, (1, 0, 2)))
            scaled = pygame.transform.smoothscale(
                source, (self.width, self.height))
            self.screen.blit(scaled, (0, 0))
        else:
            self._text("Initialising simulation...", 30, 30, self.big,
                       (220, 230, 255))
            if self.thread.error() is not None:
                self._text(str(self.thread.error()), 30, 70, self.font,
                           (255, 110, 110))

        self._draw_panel(snap)
        self._draw_hint()
        if self.show_help:
            self._draw_help()
        if time.time() < self.status_until:
            self._draw_toast()
        pygame.display.flip()

    def _text(self, text: str, x: float, y: float,
              font: pygame.font.Font, color: tuple[int, int, int]) -> None:
        self.screen.blit(font.render(text, True, color), (int(x), int(y)))

    def _draw_panel(self, snap: dict[str, Any] | None) -> None:
        lines: list[tuple[str, str]] = []
        lines.append(("N particles", f"{self.cfg.n_particles:,}"))
        lines.append(("Grid", f"{self.cfg.n_grid}³"))
        lines.append(("FPS", f"{self.clock.get_fps():5.1f}"))
        lines.append(("Sim steps/s", f"{self.thread.sim_rate:5.1f}"))
        lines.append(("Speed", f"{self.thread.speed:5.2f}x"))
        lines.append(("State", "PAUSED" if self.thread.is_paused() else "running"))
        if snap is not None:
            st = snap["stats"]
            if snap["scene"] == "collision":
                lines.append(("Time", f"{st['time_gyr']:.3f} Gyr"))
                lines.append(("Timestep", f"{st['dt_myr']:.2f} Myr"))
                lines.append(("Separation", f"{st['separation_kpc']:.1f} kpc"))
                lines.append(("Relative v", f"{st['relative_velocity']:.0f} km/s"))
                lines.append(("Since pericentre",
                              f"{st['since_pericenter_gyr']:.3f} Gyr"))
                lines.append(("t_dyn", f"{st['t_dyn_gyr']:.3f} Gyr"))
                lines.append(("t_df", f"{st['t_df_gyr']:.1f} Gyr"))
                lines.append(("Box / cell",
                              f"{st['box_kpc']:.0f} / {st['cell_kpc']:.2f} kpc"))
                lines.append(("Outside box", f"{int(st['outside'])}"))
            else:
                lines.append(("Scale factor a", f"{st['a']:.4f}"))
                lines.append(("Redshift z", f"{st['z']:.3f}"))
                lines.append(("Age", f"{st['age_gyr']:.3f} Gyr"))
                lines.append(("delta_rms", f"{st['delta_rms']:.3f}"))
                lines.append(("Nodes > 10", f"{100.0 * st['node_frac']:.2f} %"))
                lines.append(("Voids < -0.8", f"{100.0 * st['void_frac']:.1f} %"))
                lines.append(("Box / cell",
                              f"{st['box_mpch']:.0f} / {st['cell_mpch']:.2f} Mpc/h"))
                lines.append(("Slab", f"{self.slab:.2f} box"))
            if snap["scene"] == "collision":
                lines.append(("E drift", f"{100.0 * st['energy_drift']:+.3f} %"))
            lines.append(("Virial 2K/|W|", f"{st['virial']:.3f}"))
        if snap is not None and snap["scene"] == "collision":
            lines.append(("Colour", COLOR_MODES[self.color_mode]))
            lines.append(("Dark matter", "on" if self.show_dm else "off"))
        if self.recording:
            lines.append(("Recording", f"{self.record_index} frames"))

        panel_w = 330
        panel_h = 14 + 19 * len(lines)
        surface = pygame.Surface((panel_w, panel_h), pygame.SRCALPHA)
        surface.fill((6, 10, 20, 205))
        pygame.draw.rect(surface, (60, 80, 120, 180),
                         pygame.Rect(0, 0, panel_w, panel_h), 1)
        self.screen.blit(surface, (12, 12))
        self._text(self._scene_title(snap), 22, 18, self.big, (235, 242, 255))
        y = 48
        for label, value in lines:
            self._text(label, 22, y, self.small, (145, 168, 200))
            value_surface = self.small.render(value, True, (225, 235, 250))
            self.screen.blit(value_surface,
                             (panel_w - 18 - value_surface.get_width(), y))
            y += 19

    def _scene_title(self, snap: dict[str, Any] | None) -> str:
        if snap is None:
            return "GALACTIC N-BODY"
        if snap["scene"] == "collision":
            return "GALAXY COLLISION"
        return "COSMIC WEB"

    def _draw_hint(self) -> None:
        hint = ("SPACE pause   C colour   D dark matter   S save PNG   "
                "F 1920x1080   V record   T restart   H help")
        if self.cfg.scene == "cosmic":
            hint += "   [ ] slab"
        self._text(hint, 16, self.height - 26, self.small, (135, 155, 185))

    def _draw_toast(self) -> None:
        width = min(self.width - 40, 620)
        surface = pygame.Surface((width, 36), pygame.SRCALPHA)
        surface.fill((16, 28, 48, 220))
        pygame.draw.rect(surface, (90, 140, 210, 220),
                         pygame.Rect(0, 0, width, 36), 1)
        self.screen.blit(surface, (20, self.height - 72))
        self._text(self.status, 32, self.height - 66, self.small,
                   (215, 230, 250))

    def _draw_help(self) -> None:
        rows = [
            ("Mouse drag", "orbit the camera"),
            ("Mouse wheel", "zoom in / out"),
            ("Arrow keys", "orbit / change inclination"),
            ("SPACE", "pause or resume"),
            ("C", "cycle colour mode"),
            ("D", "show / hide dark matter (collision)"),
            ("+ / -", "simulation speed"),
            ("S", "save the current preview PNG"),
            ("F", "save a 1920x1080 PNG"),
            ("V", "start / stop frame recording"),
            ("[ / ]", "change the cosmic slab thickness"),
            ("R", "reset the camera"),
            ("T", "restart with a new seed"),
            ("H", "close this help"),
            ("ESC", "quit"),
        ]
        width = min(620, self.width - 60)
        height = 30 + 26 * len(rows)
        x = (self.width - width) // 2
        y = max(20, (self.height - height) // 2)
        surface = pygame.Surface((width, height), pygame.SRCALPHA)
        surface.fill((6, 10, 20, 235))
        pygame.draw.rect(surface, (80, 120, 180, 230),
                         pygame.Rect(0, 0, width, height), 1)
        self.screen.blit(surface, (x, y))
        self._text("CONTROLS", x + 20, y + 10, self.big, (235, 242, 255))
        for i, (key, action) in enumerate(rows):
            yy = y + 42 + 26 * i
            self._text(key, x + 24, yy, self.mono, (255, 205, 120))
            self._text(action, x + 190, yy, self.font, (205, 220, 240))
