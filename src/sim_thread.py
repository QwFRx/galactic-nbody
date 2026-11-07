"""Background simulation thread with lock-protected double buffering.

The UI thread never touches the live particle arrays.  The worker publishes a
finished snapshot into the inactive buffer and then flips an index under a
short lock, so a render can keep a stable reference while the next step is
already running.  This keeps the interactive loop responsive even for a
multi-million particle run.
"""

from __future__ import annotations

import threading
import time
from typing import Any

import config
from .integrator import CollisionSim, CosmicSim, make_simulation


class SimulationThread(threading.Thread):
    """Run a :class:`CollisionSim` or :class:`CosmicSim` in the background."""

    def __init__(self, cfg: config.SimConfig,
                 steps_per_second: float = 60.0) -> None:
        super().__init__(name="nbody-simulation", daemon=True)
        self.cfg = cfg
        self.base_rate = float(steps_per_second)
        self.speed = 1.0
        self.sim: CollisionSim | CosmicSim | None = None
        self.steps_done = 0
        self.sim_rate = 0.0
        self.version = 0

        self._stop_event = threading.Event()
        self._pause_event = threading.Event()
        self._ready = threading.Event()
        self._lock = threading.Lock()
        self._error: BaseException | None = None

        # Two snapshot slots: only one is written at a time.
        self._slots: list[dict[str, Any] | None] = [None, None]
        self._front = 0

        n = cfg.n_particles
        if n <= 500_000:
            self.snapshot_every = 1
        elif n <= 2_000_000:
            self.snapshot_every = 3
        else:
            self.snapshot_every = 8

    # -- worker ------------------------------------------------------------
    def run(self) -> None:
        """Create the simulation and advance it until :meth:`stop`."""
        try:
            self.sim = make_simulation(self.cfg)
            self._publish(self.sim)
            self._ready.set()
            last = time.perf_counter()
            accumulator = 0.0
            rate_window = last
            rate_steps = 0
            while not self._stop_event.is_set():
                if self._pause_event.is_set():
                    time.sleep(0.03)
                    last = time.perf_counter()
                    continue
                now = time.perf_counter()
                elapsed = min(now - last, 0.10)
                last = now
                accumulator += self.base_rate * self.speed * elapsed
                count = int(accumulator)
                if count <= 0:
                    time.sleep(0.002)
                    continue
                accumulator -= count
                count = min(count, 64)
                for _ in range(count):
                    if self._stop_event.is_set():
                        break
                    assert self.sim is not None
                    self.sim.step()
                    self.steps_done += 1
                    rate_steps += 1
                    if self.steps_done % self.snapshot_every == 0:
                        self._publish(self.sim)
                self._publish(self.sim)
                if (isinstance(self.sim, CosmicSim)
                        and self.sim.a >= 1.0 - 1.0e-10):
                    # The run has reached the present day; stop burning CPU.
                    self._pause_event.set()
                stamp = time.perf_counter()
                if stamp - rate_window >= 0.5:
                    self.sim_rate = rate_steps / (stamp - rate_window)
                    rate_window = stamp
                    rate_steps = 0
        except BaseException as exc:  # keep the message for the UI thread
            self._error = exc
        finally:
            self._ready.set()

    def _publish(self, sim: CollisionSim | CosmicSim) -> None:
        snapshot = sim.snapshot()
        with self._lock:
            back = 1 - self._front
            self._slots[back] = snapshot
            self._front = back
            self.version += 1

    # -- public API --------------------------------------------------------
    def wait_ready(self, timeout: float | None = None) -> bool:
        """Block until the first snapshot exists or construction failed."""
        return self._ready.wait(timeout)

    def error(self) -> BaseException | None:
        """Return the exception raised by the worker, if any."""
        return self._error

    def latest(self) -> dict[str, Any] | None:
        """Return the newest complete snapshot (immutable after publication)."""
        with self._lock:
            return self._slots[self._front]

    def stop(self) -> None:
        """Ask the worker to finish and wait briefly for it."""
        self._stop_event.set()
        if self.is_alive():
            self.join(timeout=2.0)

    def pause(self) -> None:
        self._pause_event.set()

    def resume(self) -> None:
        self._pause_event.clear()

    def toggle_pause(self) -> bool:
        """Flip the pause state and return ``True`` when now paused."""
        if self._pause_event.is_set():
            self._pause_event.clear()
            return False
        self._pause_event.set()
        return True

    def is_paused(self) -> bool:
        return self._pause_event.is_set()

    def set_speed(self, value: float) -> None:
        """Set the multiplicative speed in the range ``[0.1, 32]``."""
        self.speed = float(min(32.0, max(0.1, value)))

    def finished(self) -> bool:
        """True when a cosmic run has reached a = 1 (today)."""
        snap = self.latest()
        return bool(snap is not None and snap.get("scene") == "cosmic"
                    and snap.get("a", 0.0) >= 1.0 - 1.0e-9)
