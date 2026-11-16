"""Content-addressed SQLite disk cache (`cache.py`).

Guarantees deterministic replay and zero redundant API quota burn by keying on:
`sha256(endpoint + sequence_or_interval_hash + condition_hash + requested_outputs + model_version)`
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import time
from pathlib import Path
from typing import Any


class ContentAddressedCache:
    def __init__(self, cache_dir: Path | str = ".ag2_cache", model_version: str = "alphagenome-v2.0") -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = self.cache_dir / "verifier_cache.sqlite"
        self.model_version = model_version
        self._init_db()

    def _init_db(self) -> None:
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS ag2_cache (
                    cache_key TEXT PRIMARY KEY,
                    endpoint TEXT NOT NULL,
                    model_version TEXT NOT NULL,
                    created_at REAL NOT NULL,
                    payload_json TEXT NOT NULL
                )
                """
            )
            conn.commit()

    def make_key(
        self,
        endpoint: str,
        sequence_or_interval: str,
        condition_id: str,
        requested_outputs: str,
    ) -> str:
        raw = "|".join(
            [
                endpoint,
                hashlib.sha256(sequence_or_interval.encode("utf-8")).hexdigest(),
                hashlib.sha256(condition_id.encode("utf-8")).hexdigest(),
                requested_outputs,
                self.model_version,
            ]
        )
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    def get(self, cache_key: str) -> dict[str, Any] | None:
        with sqlite3.connect(self.db_path) as conn:
            row = conn.execute(
                "SELECT payload_json FROM ag2_cache WHERE cache_key = ?",
                (cache_key,),
            ).fetchone()
        if row is None:
            return None
        return json.loads(row[0])

    def put(self, cache_key: str, endpoint: str, payload: dict[str, Any]) -> None:
        payload_str = json.dumps(payload, sort_keys=True)
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                """
                INSERT OR REPLACE INTO ag2_cache
                (cache_key, endpoint, model_version, created_at, payload_json)
                VALUES (?, ?, ?, ?, ?)
                """,
                (cache_key, endpoint, self.model_version, time.time(), payload_str),
            )
            conn.commit()
