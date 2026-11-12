"""Read-only, repeatable-read database access with bounded raw-CIF batches."""
from contextlib import contextmanager
from uuid import uuid4
import pandas as pd
import psycopg
from psycopg.rows import dict_row
from .builder import build_records
from .config import Calculation, Filters, ReflectionGrid, Selection
from .errors import InsufficientDataError


class Database:
    def __init__(self, dsn="dbname=unncdb_v2", *, statement_timeout_ms=180000):
        if not isinstance(statement_timeout_ms, int) or isinstance(statement_timeout_ms, bool) or statement_timeout_ms <= 0:
            raise ValueError("statement_timeout_ms must be a positive integer")
        self._dsn = dsn
        self.statement_timeout_ms = statement_timeout_ms

    def __repr__(self):
        return "Database(<connection details hidden>)"

    @contextmanager
    def _connection(self):
        with psycopg.connect(self._dsn, connect_timeout=10, application_name="unncdb_dataset", row_factory=dict_row) as conn:
            conn.execute("SET TRANSACTION ISOLATION LEVEL REPEATABLE READ, READ ONLY")
            conn.execute("SET LOCAL search_path TO public, pg_catalog")
            conn.execute("SELECT set_config('statement_timeout', %s, true)", (str(self.statement_timeout_ms),))
            version = conn.execute("SELECT version FROM schema_migrations WHERE version=2").fetchone()
            if version is None:
                raise RuntimeError("UNNCDB schema version 2 is required")
            yield conn

    @staticmethod
    def _where(filters):
        conditions = ["src.integrity_status='verified'", "b.parse_status='parsed'",
                      "s.extraction_status='ok'", "s.is_centrosymmetric IS NOT NULL",
                      "s.space_group_num BETWEEN 1 AND 230", "s.r_factor_gt BETWEEN 0 AND %s"]
        values = [filters.r_max]
        if filters.require_verified_symmetry:
            conditions.append("s.symmetry_status='verified'")
        if filters.reject_disorder:
            conditions.append("s.disorder_status='unknown'")
        if filters.radiation_probe:
            conditions.append("lower(s.radiation_probe)=lower(%s)")
            values.append(filters.radiation_probe)
        if filters.space_groups is not None:
            conditions.append("s.space_group_num=ANY(%s)")
            values.append(list(filters.space_groups))
        if filters.max_volume is not None:
            conditions.append("s.volume<=%s")
            values.append(filters.max_volume)
        return " AND ".join(conditions), values

    @staticmethod
    def _joins():
        return (" FROM structures s JOIN cif_blocks b ON b.id=s.source_block_id "
                "JOIN cif_sources src ON src.id=b.source_id ")

    def _audit(self, conn, filters):
        where, values = self._where(filters)
        rows = conn.execute("SELECT s.space_group_num,s.is_centrosymmetric,count(*) AS candidates" +
                            self._joins() + " WHERE " + where +
                            " GROUP BY s.space_group_num,s.is_centrosymmetric ORDER BY s.space_group_num", values).fetchall()
        result = pd.DataFrame(rows, columns=["space_group_num", "is_centrosymmetric", "candidates"])
        result.attrs["note"] = "Coarse SQL candidates, before atom checks, deduplication and final class balancing. Unknown disorder is not confirmed order."
        return result

    def audit(self, filters=None):
        """Return candidate counts per space group; final acceptance requires build()."""
        filters = filters or Filters()
        with self._connection() as conn:
            return self._audit(conn, filters)

    def _records(self, conn, filters, selection, batch_size):
        where, params = self._where(filters)
        statement = ("SELECT s.id AS structure_id,s.refcode,s.space_group_num,s.is_centrosymmetric,"
                     "s.a,s.b,s.c,s.alpha,s.beta,s.gamma,s.r_factor_gt,b.id AS source_block_id,"
                     "src.id AS source_id,src.sha256 AS source_sha256,b.sha256 AS cif_sha256,"
                     "i.importer_version,i.parser_version" + self._joins() +
                     " JOIN cif_imports i ON i.id=b.import_id WHERE " + where +
                     " ORDER BY md5(b.sha256 || ':' || %s),s.id")
        params.append(str(selection.seed))
        if selection.max_candidates is not None:
            statement += " LIMIT %s"
            params.append(selection.max_candidates)
        with conn.cursor(name="dataset_" + uuid4().hex) as cursor:
            cursor.itersize = batch_size
            cursor.execute(statement, params)
            while batch := cursor.fetchmany(batch_size):
                # Keep TOASTed CIF bytes out of the global sort; fetch only one batch.
                ids = [row["source_block_id"] for row in batch]
                raw = conn.execute("SELECT id,raw_cif FROM cif_blocks WHERE id=ANY(%s)", (ids,)).fetchall()
                blobs = {row["id"]: bytes(row["raw_cif"]) for row in raw}
                for record in batch:
                    yield record, blobs[record["source_block_id"]]

    def build(self, n=1000, *, filters=None, grid=None, calculation=None, balance=True, seed=42,
              one_per_family=True, max_per_space_group=None, max_candidates=None,
              allow_partial=False, batch_size=128, progress=None):
        """Create a Dataset in memory. No database writes or train/test split.

        A failed exact request raises InsufficientDataError with a report and
        partial_dataset; allow_partial=True returns an explicitly incomplete result.
        """
        if not isinstance(batch_size, int) or isinstance(batch_size, bool) or batch_size <= 0:
            raise ValueError("batch_size must be a positive integer")
        if progress is not None and not callable(progress):
            raise ValueError("progress must be a callable or None")
        filters, grid, calculation = filters or Filters(), grid or ReflectionGrid.cube(4), calculation or Calculation()
        selection = Selection(n=n, balance=balance, seed=seed, one_per_family=one_per_family,
                              max_per_space_group=max_per_space_group, max_candidates=max_candidates,
                              allow_partial=allow_partial)
        with self._connection() as conn:
            audit = self._audit(conn, filters)
            counts = {0: 0, 1: 0}
            for row in audit.itertuples(index=False):
                counts[int(row.is_centrosymmetric)] += int(row.candidates)
            capacity = dict(counts)
            if max_per_space_group is not None:
                capacity = {0: 0, 1: 0}
                for row in audit.itertuples(index=False):
                    capacity[int(row.is_centrosymmetric)] += min(int(row.candidates), max_per_space_group)
            if not allow_partial and ((balance and min(capacity.values()) < n // 2) or sum(capacity.values()) < n):
                raise InsufficientDataError(f"Too few candidates under SQL/group-cap conditions: {capacity}",
                                            {"coarse_class_counts": counts,
                                             "class_capacity": capacity, "requested": n})
            tx = conn.execute("SELECT txid_current_snapshot()::text AS snapshot,current_database() AS name").fetchone()
            provenance = {"name": tx["name"], "snapshot": tx["snapshot"], "schema_version": 2,
                          "isolation": "repeatable read, read only", "coarse_class_counts": counts,
                          "class_capacity": capacity,
                          "candidate_order": "md5(cif_sha256 + ':' + seed), structure_id"}
            return build_records(self._records(conn, filters, selection, batch_size),
                                 selection=selection, filters=filters, grid=grid,
                                 calculation=calculation, provenance=provenance, progress=progress)
