"""
UniProt Reference Set Retrieval Library
========================================
Retrieves protein sequences from the local UniProt database.

Can be used in two ways:

1. As a command-line tool (for quick retrieval without writing custom scripts):
   python get_reference_uniprot_set_lib.py -version 2026_01 -taxonomy 9606 10090

2. As an importable library in your scripts:

   # --- Quickstart (one-liner) ---
   from get_reference_uniprot_set_lib import fetch_sequences

   records = fetch_sequences(version="2026_01", taxon_ids=[9606, 10090])
   for r in records:
       print(r["accession"], r["organism"])

   # --- HMM hit retrieval (one-liner) ---
   from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit

   records = fetch_sequences_by_hmm_hit(
       version="2026_01", hmm_query="Homeodomain", taxon_ids=[9606]
   )

   # --- Context manager (recommended for multiple queries) ---
   from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config

   with UniProtRetriever(get_db_config()) as db:
       records    = db.get_proteins(version="2026_01", taxon_ids=[9606, 10090])
       hmm_recs   = db.get_proteins_by_hmm_hit(version="2026_01", hmm_query="Homeodomain")
       fasta_str  = db.to_fasta_string(records)       # in-memory FASTA string
       seqrecords = db.to_biopython(records)           # list of BioPython SeqRecord
       db.export_fasta(records, "2026_01", "mammals", output_dir="./fastas")
"""

import mysql.connector
import argparse
import os
import sys
from dotenv import load_dotenv
import getpass
import pandas as pd

load_dotenv()

# Public API — what users get when they do `from ... import *`
__all__ = [
    "UniProtRetriever",
    "get_db_config",
    "fetch_sequences",
    "fetch_sequences_by_hmm_hit",
    "fetch_fasta_string",
    "fetch_fasta_string_by_hmm_hit",
    "fetch_domains_by_accession",
    "fetch_sequences_by_accession",
    "fetch_domains_by_go",
    "fetch_presence_absence_matrix",
    "fetch_accessions_for_cell",
    "fetch_subprofile_hits",
    "fetch_domain_architectures",
    "fetch_accessions_with_taxids",
    "fetch_highres_profile",
]


# ==========================================================================================================

#                                       CONFIGURATION HELPER

# ==========================================================================================================


def get_db_config(host=None, user=None, password=None, database=None):
    """
    Build the database config dict, falling back to environment variables
    and then to the lab defaults.

    users can override any value when calling this function, or by
    setting environment variables in a .env file:
        DB_HOST, DB_USER, DB_PASSWORD, DB_NAME
    """
    return {
        "host": host or os.getenv("DB_HOST", "localhost"),
        "user": user or os.getenv("DB_USER", "user"),
        "password": password or os.getenv("DB_PASSWORD", "your_password"),
        "database": database or os.getenv("DB_NAME", "uniprot_db"),
    }


# ==========================================================================================================

#                                               CORE CLASS
# ==========================================================================================================


class UniProtRetriever:
    """
    Retrieves UniProt reference sets from the local CGLab database.

    Supports filtering by TaxID, Proteome ID, Pfam domain, GO term, and HMM search
    results from the local Pfam-A hmmsearch conducted by pyhmmer hmmsearch.

    Results can be returned as raw dicts, a FASTA string, BioPython
    SeqRecord objects, or written directly to a .fasta file.

    Recommended usage — context manager (handles connect/close):

        with UniProtRetriever(get_db_config()) as db:
            records = db.get_proteins(version="2026_01", taxon_ids=[9606])
            seqs = db.to_biopython(records)

        # HMM-based retrieval for tree building:
        with UniProtRetriever(get_db_config()) as db:
            records = db.get_proteins_by_hmm_hit(
                version="2026_01", hmm_query="Homeodomain", taxon_ids=[9606]
            )
            seqs = db.to_biopython(records)
    """

    def __init__(self, config):

        self.config = config
        self.conn = None
        self.cursor = None

    # ==========================================================================================================

    #                                       CONNECTION MANAGEMENT

    # ==========================================================================================================
    def connect(self):
        """
        Open the database connection.

        """
        self.conn = mysql.connector.connect(**self.config)
        self.cursor = self.conn.cursor(dictionary=True)

    def close(self):
        """Close the database connection (safe to call even if not connected)."""
        if self.conn and self.conn.is_connected():
            self.cursor.close()
            self.conn.close()

    # Context-manager support — in order to be able to use `with UniProtRetriever(...) as db:`
    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    # ==========================================================================================================

    #                                       UTILITY QUERIES

    # ==========================================================================================================

    def list_available_versions(self):
        """
        Print a summary table of all versions present in the database.

        """
        query = """
            SELECT version,
                   COUNT(*)                  AS protein_count,
                   COUNT(DISTINCT taxon_id)  AS taxon_count,
                   COUNT(DISTINCT proteome_id) AS proteome_count
            FROM   proteins
            GROUP  BY version
            ORDER  BY version DESC
        """
        self.cursor.execute(query)
        versions = self.cursor.fetchall()
        if versions:
            print("\nAvailable versions in database:")
            print("-" * 70)
            for v in versions:
                print(
                    f"Version: {v['version']} | "
                    f"Proteins: {v['protein_count']:,} | "
                    f"Taxa: {v['taxon_count']:,} | "
                    f"Proteomes: {v['proteome_count']:,}"
                )
            print("-" * 70)
        else:
            print("\nNo versions found in database.")
        return versions

    def get_proteome_ids(self, version):
        """
        Return all unique Proteome IDs present for a given version.

        Parameters
        ----------
        version : str  e.g. "2026_01"

        Returns
        -------
        list[str]
        """
        query = "SELECT DISTINCT proteome_id FROM proteins WHERE version = %s"
        self.cursor.execute(query, (version,))
        return [row["proteome_id"] for row in self.cursor.fetchall()]

    # ==========================================================================================================

    #                                               MAIN RETRIEVAL

    # ==========================================================================================================

    def get_proteins(
        self,
        version,
        taxon_ids=None,
        proteome_id=None,
        go_id=None,
        pfam_id=None,
    ):
        """
        Retrieve protein records from the database.

        All filters are optional and combinable.


        """
        query = """
            SELECT p.accession, p.name, p.organism,
                   p.taxon_id, p.proteome_id, s.sequence
            FROM   proteins  p
            JOIN   sequences s ON p.seq_id = s.seq_id
        """
        if go_id:
            query += (
                " JOIN protein_go pg"
                " ON p.accession = pg.accession AND p.version = pg.version"
            )
        if pfam_id:
            query += (
                " JOIN protein_pfam pp"
                " ON p.accession = pp.accession AND p.version = pp.version"
            )

        where_clauses = ["p.version = %s"]
        params = [version]

        if taxon_ids is not None:
            if isinstance(taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(taxon_ids))
                where_clauses.append(f"p.taxon_id IN ({placeholders})")
                params.extend(taxon_ids)
            else:
                where_clauses.append("p.taxon_id = %s")
                params.append(taxon_ids)

        if proteome_id:
            where_clauses.append("p.proteome_id = %s")
            params.append(proteome_id)

        if go_id:
            where_clauses.append("pg.go_id = %s")
            params.append(go_id)

        if pfam_id:
            where_clauses.append("pp.pfam_id = %s")
            params.append(pfam_id)

        query += " WHERE " + " AND ".join(where_clauses)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ==========================================================================================================

    #                                       HMM HIT BASED RETRIEVAL

    # ==========================================================================================================
    def get_proteins_by_hmm_hit(
        self,
        version,
        hmm_query,
        evalue_cutoff=None,
        taxon_ids=None,
        exclude_taxon_ids=None,
    ):
        """
        Fetch sequences for all proteins with a hit to a given HMM profile.
        Suitable for phylogenetic tree building pipelines.

        Queries the local hmm_search_results table that contains results from a full
        Pfam-A hmmsearch using gathering thresholds. Results are already
        filtered at the profile-specific trusted cutoff. So, evalue_cutoff
        provides an additional filter on top of that.

        Results are deduplicated: a protein with five Homeodomain repeats
        appears only once in the output, since only its sequence is needed.

        Suitable for phylogenetic tree building pipelines.


        Examples
        --------
        # All Homeodomain proteins across all taxa:
        records = db.get_proteins_by_hmm_hit("2026_01", "Homeodomain")

        # Human + mouse kinases with strict E-value:
        records = db.get_proteins_by_hmm_hit(
            "2026_01", "PF00069", evalue_cutoff=1e-10, taxon_ids=[9606, 10090]
        )
        seqrecords = db.to_biopython(records)
        db.export_fasta(records, "2026_01", "kinase_human_mouse", output_dir="./fastas")
        """

        # Uses DISTINCT because a protein might have 5 repeating Homeodomains, but we only need to fetch its sequence once.
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""
        query = f"""
            SELECT DISTINCT p.accession, p.name, p.organism,
                            p.taxon_id, p.proteome_id, s.sequence
            FROM   hmm_search_results h
            JOIN   proteins p  ON h.accession = p.accession AND h.version = p.version
            JOIN   sequences s ON p.seq_id = s.seq_id
            WHERE  h.version = %s
              {evalue_clause}
              AND  (h.hmm_name = %s OR h.hmm_accession LIKE %s)
        """

        # We use LIKE %s with a wildcard (%) so if the user searches for "PF00046",
        # it successfully matches "PF00046.36" in your database.
        params = [version]
        if evalue_cutoff is not None:
            params.append(evalue_cutoff)
        params += [hmm_query, f"{hmm_query}%"]

        if taxon_ids is not None:
            if isinstance(taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(taxon_ids))
                query += f" AND h.taxon_id IN ({placeholders})"
                params.extend(taxon_ids)
            else:
                query += " AND h.taxon_id = %s"
                params.append(taxon_ids)

        if exclude_taxon_ids is not None:
            if isinstance(exclude_taxon_ids, (list, tuple)):
                placeholders = ", ".join(["%s"] * len(exclude_taxon_ids))
                query += f" AND h.taxon_id NOT IN ({placeholders})"
                params.extend(exclude_taxon_ids)

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ==========================================================================================================

    #                                       DOMAIN LOOKUP BY ACCESSION

    # ==========================================================================================================

    def get_domains_by_accession(
        self,
        version,
        accessions,
        evalue_cutoff=None,
        chunk_size=5000,
    ):
        """
        Return all Pfam domain hits from hmm_search_results for one or
        more protein accessions.

        One row is returned per domain occurrence — a protein with five
        Ankyrin repeats produces five rows, one per repeat. This is
        intentional: use this method when you need domain coordinates and
        scores, not just the sequence.

        Large accession lists are chunked automatically (default 5,000
        per query) so the IN (...) clause never becomes a bottleneck.


        Examples
        --------
        # Single protein — all its Pfam domains
        domains = db.get_domains_by_accession("2026_01", "P04637")

        # Multiple proteins
        domains = db.get_domains_by_accession(
            "2026_01", ["P04637", "P10275", "Q8BFR5"]
        )
        for d in domains:
            print(f"{d['accession']}  {d['hmm_name']}  "
                  f"{d['ali_from']}–{d['ali_to']}  E={d['full_evalue']:.2e}")
        """
        if isinstance(accessions, str):
            accessions = [accessions]
        if not accessions:
            return []

        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""

        all_results = []
        try:
            for i in range(0, len(accessions), chunk_size):
                chunk = accessions[i : i + chunk_size]
                placeholders = ", ".join(["%s"] * len(chunk))

                query = f"""
                    SELECT
                        h.accession,
                        p.name          AS protein_name,
                        p.organism,
                        p.taxon_id,
                        p.proteome_id,
                        h.hmm_name,
                        h.hmm_accession,
                        h.domain_number,
                        h.domain_count,
                        h.full_evalue,
                        h.full_score,
                        h.domain_evalue,
                        h.domain_score,
                        h.ali_from,
                        h.ali_to,
                        h.hmm_from,
                        h.hmm_to,
                        h.env_from,
                        h.env_to
                    FROM  hmm_search_results h
                    JOIN  proteins p
                    ON  h.accession = p.accession
                    AND  h.version   = p.version
                    WHERE h.version = %s
                    {'AND h.full_evalue <= %s' if evalue_cutoff is not None else ''}
                    AND (h.accession IN ({placeholders}) OR p.name IN ({placeholders}))
                    ORDER BY h.accession, h.ali_from
                """

                params = [version]
                if evalue_cutoff is not None:
                    params.append(evalue_cutoff)
                params += chunk + chunk

                self.cursor.execute(query, tuple(params))
                all_results.extend(self.cursor.fetchall())
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        return all_results

    def get_domains_by_go(self, version, go_id, evalue_cutoff=None):
        """
        Given a GO term, return all HMM profiles found in proteins annotated
        with that term, with counts.
        """
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""
        query = f"""
            SELECT h.hmm_name, h.hmm_accession,
                COUNT(DISTINCT p.accession) AS protein_count
            FROM   protein_go pg
            JOIN   proteins p ON pg.accession = p.accession AND pg.version = p.version
            JOIN   hmm_search_results h ON p.accession = h.accession AND p.version = h.version
            WHERE  pg.version = %s
            AND    pg.go_id = %s
            {evalue_clause}
            GROUP  BY h.hmm_name, h.hmm_accession
            ORDER  BY protein_count DESC
        """
        params = [version, go_id]
        if evalue_cutoff is not None:
            params.append(evalue_cutoff)
        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ------------------------------------------------------------------
    # Presence / Absence matrix  (Step 1)
    # ------------------------------------------------------------------

    def get_presence_absence_matrix(
        self,
        version,
        pfam_queries,
        taxon_ids=None,
        evalue_cutoff=None,
    ):
        """
        Build a presence/absence matrix: for each (taxon, Pfam profile) pair,
        count how many distinct proteins carry that profile.

        Step 1 of the two-tier resolution workflow.
        Returns flat rows — pivot with pandas in the UI layer.


        """
        if isinstance(pfam_queries, str):
            pfam_queries = [pfam_queries]
        if not pfam_queries:
            return []

        pfam_conditions = " OR ".join(
            ["(h.hmm_name = %s OR h.hmm_accession LIKE %s)"] * len(pfam_queries)
        )
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""

        # ── Query 1: count distinct proteins per (taxon, profile) ────────
        # No JOIN to proteomes here — just hmm_search_results.
        query = f"""
            SELECT
                h.taxon_id,
                h.hmm_name,
                h.hmm_accession,
                h.hmm_type,
                COUNT(DISTINCT h.accession) AS protein_count
            FROM   hmm_search_results h
            WHERE  h.version = %s
              {evalue_clause}
              AND ({pfam_conditions})
        """

        params = [version]
        if evalue_cutoff is not None:
            params.append(evalue_cutoff)
        for q in pfam_queries:
            params.extend([q, f"{q}%"])

        if taxon_ids is not None:
            if isinstance(taxon_ids, (int, str)):
                taxon_ids = [taxon_ids]
            placeholders = ", ".join(["%s"] * len(taxon_ids))
            query += f" AND h.taxon_id IN ({placeholders})"
            params.extend(taxon_ids)

        query += """
            GROUP BY h.taxon_id, h.hmm_name, h.hmm_accession, h.hmm_type
            ORDER BY h.taxon_id, protein_count DESC
        """

        try:
            self.cursor.execute(query, tuple(params))
            rows = self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        if not rows:
            return []

        # ── Query 2: resolve taxon_id → organism name ────────────────────
        # We collect the unique taxon_ids from the results, then ask the
        # proteins table for one representative organism string per taxon.
        # MIN(organism) is arbitrary but deterministic — every protein from
        # the same taxon has the same organism string anyway.
        found_taxon_ids = list({r["taxon_id"] for r in rows})
        placeholders = ", ".join(["%s"] * len(found_taxon_ids))

        name_query = f"""
            SELECT   taxon_id, MIN(organism) AS scientific_name
            FROM     proteins
            WHERE    version = %s
              AND    taxon_id IN ({placeholders})
            GROUP BY taxon_id
        """
        try:
            self.cursor.execute(name_query, tuple([version] + found_taxon_ids))
            name_rows = self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(
                f"Database query failed (taxon name lookup): {err}"
            ) from err

        taxon_name_map = {r["taxon_id"]: r["scientific_name"] for r in name_rows}

        # ── Attach scientific_name to each result row ────────────────────
        for row in rows:
            row["scientific_name"] = taxon_name_map.get(row["taxon_id"])

        return rows

    # ------------------------------------------------------------------
    # Accession retrieval for one matrix cell  (Step 1 → Step 2 bridge)
    # ------------------------------------------------------------------

    def get_accessions_for_cell(
        self,
        version,
        pfam_query,
        taxon_id,
        evalue_cutoff=None,
    ):
        """
        Return the accessions (and basic metadata) for a single cell
        of the presence/absence matrix.

        Call this when the user clicks a cell in the UI. The returned
        accessions are then passed directly into get_subprofile_hits()
        and/or get_domain_architectures() for step-2 drill-down.

        Example
        -------
        # User clicks the (9606, "Homeodomain") cell:
        accs = db.get_accessions_for_cell("2026_01", "Homeodomain", 9606)
        accession_list = [r["accession"] for r in accs]
        # → pass accession_list into get_subprofile_hits() or get_domain_architectures()
        """
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""

        query = f"""
            SELECT DISTINCT
                h.accession,
                h.protein_name,
                p.organism,
                h.taxon_id,
                h.proteome_id
            FROM   hmm_search_results h
            JOIN   proteins p ON h.accession = p.accession AND h.version = p.version
            WHERE  h.version = %s
              {evalue_clause}
              AND  (h.hmm_name = %s OR h.hmm_accession LIKE %s)
              AND  h.taxon_id = %s
            ORDER BY h.accession
        """

        params = [version]
        if evalue_cutoff is not None:
            params.append(evalue_cutoff)
        params.extend([pfam_query, f"{pfam_query}%", int(taxon_id)])

        try:
            self.cursor.execute(query, tuple(params))
            return self.cursor.fetchall()
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

    # ------------------------------------------------------------------
    # Sub-profile hits  (Step 2, Path A — deeper HMM resolution)
    # ------------------------------------------------------------------

    def get_subprofile_hits(
        self,
        version,
        accessions,
        evalue_cutoff=None,
        exclude_queries=None,
        chunk_size=5000,
    ):
        """
        Given a set of accessions (one matrix cell), return ALL HMM profiles
        that hit at least one of those proteins — sorted by how many proteins
        carry each profile.

        Step 2, Path A: "deeper HMM resolution."
        The original broad Pfam (e.g. Homeodomain) will appear at the top.
        Below it you will find more specific profiles — TIGRFAMs, sub-Pfams,
        SUPERFAMILYs — that subdivide this protein set into functional
        sub-groups (e.g. Histamine DC, Tyramine DC, Dopamine DC).

        Returns:

            Sorted by protein_count descending. Keys:
                hmm_name        str
                hmm_accession   str    (versioned, e.g. "TIGR00001.1")
                hmm_type        str    "Pfam" | "TIGRFAM" | "SUPERFAMILY" | "other"
                protein_count   int    how many of the input accessions carry it
                coverage        float  protein_count / total_input_accessions
                best_evalue     float  lowest (best) full_evalue in the set
                best_score      float  highest full_score in the set

        Example
        -------
        accs = [r["accession"] for r in
                db.get_accessions_for_cell("2026_01", "Homeodomain", 9606)]

        hits = db.get_subprofile_hits(
            "2026_01", accs,
            exclude_queries=["Homeodomain", "PF00046"],
        )
        # hits[0] → most common co-occurring profile in this human Homeodomain set
        """
        if isinstance(accessions, str):
            accessions = [accessions]
        if not accessions:
            return []

        total = len(accessions)
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""

        all_rows = []
        try:
            for i in range(0, len(accessions), chunk_size):
                chunk = accessions[i : i + chunk_size]
                placeholders = ", ".join(["%s"] * len(chunk))

                query = f"""
                    SELECT
                        h.hmm_name,
                        h.hmm_accession,
                        h.hmm_type,
                        COUNT(DISTINCT h.accession) AS protein_count,
                        MIN(h.full_evalue)          AS best_evalue,
                        MAX(h.full_score)           AS best_score
                    FROM   hmm_search_results h
                    WHERE  h.version = %s
                      {evalue_clause}
                      AND  h.accession IN ({placeholders})
                    GROUP BY h.hmm_name, h.hmm_accession, h.hmm_type
                """

                params = [version]
                if evalue_cutoff is not None:
                    params.append(evalue_cutoff)
                params.extend(chunk)

                self.cursor.execute(query, tuple(params))
                all_rows.extend(self.cursor.fetchall())

        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        # Merge chunks: a profile might appear in multiple chunks,
        # so aggregate across chunks in Python.
        merged = {}
        for row in all_rows:
            key = (row["hmm_name"], row["hmm_accession"])
            if key not in merged:
                merged[key] = dict(row)
            else:
                merged[key]["protein_count"] += row["protein_count"]
                merged[key]["best_evalue"] = min(
                    merged[key]["best_evalue"], row["best_evalue"]
                )
                merged[key]["best_score"] = max(
                    merged[key]["best_score"], row["best_score"]
                )

        results = list(merged.values())

        # Add coverage fraction: what fraction of the input set carries each profile?
        for r in results:
            r["coverage"] = round(r["protein_count"] / total, 4)

        # Optionally hide the original broad query from the output so the
        # sub-profiles are immediately visible without scrolling past the hit
        # that was already shown in step 1.
        if exclude_queries:
            exclude_set = set(exclude_queries)
            results = [
                r
                for r in results
                if r["hmm_name"] not in exclude_set
                and not any(r["hmm_accession"].startswith(q) for q in exclude_set)
            ]

        results.sort(key=lambda r: r["protein_count"], reverse=True)
        return results

    # ------------------------------------------------------------------
    # Domain architecture patterns  (Step 2, Path B — domain co-occurrence)
    # ------------------------------------------------------------------

    def get_domain_architectures(
        self,
        version,
        accessions,
        evalue_cutoff=None,
        collapse_repeats=True,
        chunk_size=5000,
    ):
        """
        For a set of accessions, characterise each protein's full domain
        architecture and count how many proteins share each pattern.

        This is Step 2, Path B: "domain co-occurrence resolution."
        Instead of asking which sub-HMMs are enriched (Path A), this asks:
        what combinations of domains define this protein set?

        Architecture is defined as the ordered, left-to-right sequence of
        domain names on a protein (ordered by ali_from, i.e. alignment start
        in the protein sequence).

        Example
        -------
        accs = [r["accession"] for r in
                db.get_accessions_for_cell("2026_01", "Homeodomain", 9606)]

        archs = db.get_domain_architectures("2026_01", accs)
        for a in archs[:5]:
            print(a["architecture"], "→", a["protein_count"], "proteins")
        # Homeodomain+PBC → 38 proteins
        # Homeodomain     → 12 proteins
        # TALE+Homeodomain → 4 proteins
        """
        if isinstance(accessions, str):
            accessions = [accessions]
        if not accessions:
            return []

        total = len(accessions)
        evalue_clause = "AND h.full_evalue <= %s" if evalue_cutoff is not None else ""

        # Fetch all domain hits for these accessions, ordered by position.
        # One row per domain occurrence — a protein with 5 ANK repeats gives 5 rows.
        # We collect them in Python and build architecture strings per accession.
        all_rows = []
        try:
            for i in range(0, len(accessions), chunk_size):
                chunk = accessions[i : i + chunk_size]
                placeholders = ", ".join(["%s"] * len(chunk))

                query = f"""
                    SELECT
                        h.accession,
                        h.hmm_name,
                        h.hmm_accession,
                        h.ali_from,
                        h.ali_to
                    FROM   hmm_search_results h
                    WHERE  h.version = %s
                      {evalue_clause}
                      AND  h.accession IN ({placeholders})
                    ORDER BY h.accession, h.ali_from
                """

                params = [version]
                if evalue_cutoff is not None:
                    params.append(evalue_cutoff)
                params.extend(chunk)

                self.cursor.execute(query, tuple(params))
                all_rows.extend(self.cursor.fetchall())

        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        # Group domain rows by accession, preserving ali_from order (already sorted).
        # Build architecture string and arch_accessions string per protein.
        from collections import defaultdict

        domains_per_acc = defaultdict(list)  # accession → [(hmm_name, hmm_accession)]
        for row in all_rows:
            domains_per_acc[row["accession"]].append(
                (row["hmm_name"], row["hmm_accession"])
            )

        # Proteins with NO hits at all (passed evalue filter or simply no domains)
        # will be absent from domains_per_acc — they get architecture "—".
        arch_counts = defaultdict(
            lambda: {"protein_count": 0, "examples": [], "arch_accessions": ""}
        )

        for acc in accessions:
            domain_list = domains_per_acc.get(acc, [])

            if not domain_list:
                arch_str = "—"
                arch_acc_str = "—"
            else:
                if collapse_repeats:
                    # Keep only first occurrence of each domain name,
                    # preserving positional order.
                    seen_names = []
                    seen_accs = []
                    seen_set = set()
                    for name, hacc in domain_list:
                        if name not in seen_set:
                            seen_set.add(name)
                            seen_names.append(name)
                            seen_accs.append(hacc)
                else:
                    seen_names = [d[0] for d in domain_list]
                    seen_accs = [d[1] for d in domain_list]

                arch_str = "+".join(seen_names)
                arch_acc_str = "+".join(seen_accs)

            entry = arch_counts[arch_str]
            entry["protein_count"] += 1
            entry["arch_accessions"] = (
                arch_acc_str  # same for all proteins sharing arch
            )
            if len(entry["examples"]) < 5:
                entry["examples"].append(acc)

        # Format output
        results = []
        for arch_str, data in arch_counts.items():
            results.append(
                {
                    "architecture": arch_str,
                    "arch_accessions": data["arch_accessions"],
                    "protein_count": data["protein_count"],
                    "coverage": round(data["protein_count"] / total, 4),
                    "example_accessions": data["examples"],
                }
            )

        results.sort(key=lambda r: r["protein_count"], reverse=True)
        return results

    def get_proteins_by_accession(
        self,
        version,
        accessions,
        chunk_size=5000,
    ):
        """
        Fetch full protein records (including sequence) for one or more
        UniProt accessions.

        Large accession lists are chunked automatically (default 5,000
        per query) so the IN (...) clause never happens a bottleneck.


        Examples
        --------
        # Single protein
        records = db.get_proteins_by_accession("2026_01", "P04637")

        # Multiple proteins
        records = db.get_proteins_by_accession(
            "2026_01", ["P04637", "P10275", "Q8BFR5"]
        )
        fasta_str = db.to_fasta_string(records)
        """
        if isinstance(accessions, str):
            accessions = [accessions]
        if not accessions:
            return []

        all_results = []
        try:
            for i in range(0, len(accessions), chunk_size):
                chunk = accessions[i : i + chunk_size]
                placeholders = ", ".join(["%s"] * len(chunk))

                query = f"""
                    SELECT p.accession, p.name, p.organism,
                           p.taxon_id, p.proteome_id, s.sequence
                    FROM   proteins  p
                    JOIN   sequences s ON p.seq_id = s.seq_id
                    WHERE  p.version = %s
                      AND  (p.accession IN ({placeholders}) OR p.name IN ({placeholders}))
                """

                params = [version] + chunk + chunk
                self.cursor.execute(query, tuple(params))
                all_results.extend(self.cursor.fetchall())
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        return all_results

    def get_accessions_with_taxids(
        self,
        version,
        accessions,
        chunk_size=5000,
    ):
        """Bulk lookup: accession → taxon_id + scientific name.

        Used to map tree leaves (UniProt accessions) back to the taxa
        that own them, so subclade memberships can be turned into the
        rows of a high-resolution phylogenetic profile matrix.

        Large accession lists are chunked automatically (default 5,000
        per query) so the IN (...) clause never becomes a bottleneck.


        Examples
        --------
        rows = db.get_accessions_with_taxids("2026_01", ["P04637", "P10275"])
        acc_to_taxon = {r["accession"]: r["taxon_id"] for r in rows}
        """

        if isinstance(accessions, str):
            accessions = [accessions]
        if not accessions:
            return []

        # Deduplicate — trees can have repeated leaves after some ops,
        # so the same accession can show up across subclades by mistake.
        accessions = list(set(accessions))

        all_results = []
        try:
            for i in range(0, len(accessions), chunk_size):
                chunk = accessions[i : i + chunk_size]
                placeholders = ", ".join(["%s"] * len(chunk))

                query = f"""
                    SELECT   accession,
                             taxon_id,
                             organism AS scientific_name
                    FROM     proteins
                    WHERE    version = %s
                      AND    accession IN ({placeholders})
                """
                params = [version] + chunk
                self.cursor.execute(query, tuple(params))
                all_results.extend(self.cursor.fetchall())
        except mysql.connector.Error as err:
            raise RuntimeError(f"Database query failed: {err}") from err

        return all_results

    # High-resolution phylogenetic profile (gets the output of the subclade partition)
    def get_highres_profile(
        self,
        version,
        pfam_subclade_map,
        taxon_ids=None,
        binary=False,
    ):
        """
        Build a high-resolution phylogenetic profile matrix.

        Per-Pfam gene-tree subclades become the columns of the matrix.
        Each cell counts how many proteins of that Pfam-subclade are
        carried by that taxon. With `binary=True`, cells become 0/1
        instead of counts.

         - Column order is preserved exactly as in the input
            Don't rely on alphabetical order.
        - `missing_accessions` is not an error — common causes are
          (a) tree was built against a different UniProt version, or
          (b) the accession was renamed/removed between versions.

        """
        # 1. Flatten input --> set of every unique accession
        all_accessions = set()
        for pfam, subs in pfam_subclade_map.items():
            for label, accs in subs.items():
                all_accessions |= set(accs)

        empty_return = {
            "matrix": pd.DataFrame(),
            "taxon_names": {},
            "missing_accessions": set(),
            "column_origin": {},
        }

        if not all_accessions:
            return empty_return

        # 2. Bulk taxid lookup
        rows = self.get_accessions_with_taxids(version, list(all_accessions))
        acc_to_taxon = {r["accession"]: r["taxon_id"] for r in rows}
        taxon_names = {r["taxon_id"]: r["scientific_name"] for r in rows}

        missing = all_accessions - set(acc_to_taxon.keys())

        # 3. Build long-form (taxon, column) records
        column_origin = {}
        long_records = []
        for pfam, subs in pfam_subclade_map.items():
            for label, accs in subs.items():
                col = f"{pfam}-{label}"
                column_origin[col] = (pfam, label)
                for acc in accs:
                    tx = acc_to_taxon.get(acc)
                    if tx is None:
                        continue  # accession not in DB --> already tracked in `missing`
                    long_records.append((tx, col))

        if not long_records:
            empty_return["taxon_names"] = taxon_names
            empty_return["missing_accessions"] = missing
            empty_return["column_origin"] = column_origin
            return empty_return

        long_df = pd.DataFrame(long_records, columns=["taxon_id", "column"])

        # 4. Pivot --> wide matrix, counts
        matrix = long_df.groupby(["taxon_id", "column"]).size().unstack(fill_value=0)

        # 5. All columns present (even empty) and preserve order
        for col in column_origin:
            if col not in matrix.columns:
                matrix[col] = 0
        matrix = matrix[list(column_origin.keys())]

        # 6. Optional row filter / reorder
        if taxon_ids is not None:
            if isinstance(taxon_ids, (int, str)):
                taxon_ids = [taxon_ids]
            taxon_ids = [int(t) for t in taxon_ids]
            matrix = matrix.reindex(taxon_ids, fill_value=0)

        if binary:
            matrix = (matrix > 0).astype(int)

        # Name the axes for cleaner downstream display
        matrix.index.name = "taxon_id"
        matrix.columns.name = "pfam_subclade"

        return {
            "matrix": matrix,
            "taxon_names": taxon_names,
            "missing_accessions": missing,
            "column_origin": column_origin,
        }

    # ------------------------------------------------------------------
    # Output helpers
    # ------------------------------------------------------------------

    def to_fasta_string(self, records):
        """
        Convert a list of protein records to a FASTA-formatted string.

        This does NOT write any file — the string is returned so callers
        can pass it directly to parsers (e.g. Bio.SeqIO.parse via StringIO).

        """
        lines = []
        for rec in records:
            header = f">{rec['taxon_id']}.{rec['accession']}"
            lines.append(f"{header}\n{rec['sequence']}")
        return "\n".join(lines) + "\n"

    def to_biopython(self, records):
        """
        Convert records to a list of BioPython SeqRecord objects.

        Requires BioPython (pip install biopython).
        Useful for alignment tools, tree building (ETE3, DendroPy, etc.),
        and any workflow that expects SeqRecord objects.

        """
        try:
            from Bio.Seq import Seq
            from Bio.SeqRecord import SeqRecord
        except ImportError as e:
            raise ImportError(
                "BioPython is required for to_biopython(). "
                "Install it with:  pip install biopython"
            ) from e

        seqrecords = []
        for rec in records:
            sr = SeqRecord(
                Seq(rec["sequence"]),
                id=f"{rec['taxon_id']}.{rec['accession']}",
                name=rec["accession"],
                description=f"{rec['name']} [{rec['organism']}] "
                f"UP={rec['proteome_id']}",
            )
            seqrecords.append(sr)
        return seqrecords

    def export_fasta(
        self, records, version, identifier, output_dir=None, filename=None
    ):
        """
        Write protein records to a FASTA file.

        """
        if filename is None:
            filename = f"uniprot_{identifier}_{version}.fasta"

        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            filepath = os.path.join(output_dir, filename)
        else:
            filepath = filename

        with open(filepath, "w") as f:
            f.write(self.to_fasta_string(records))

        print(f" Exported {len(records):,} sequences → {os.path.abspath(filepath)}")
        return os.path.abspath(filepath)


# ==========================================================================================================

#                            MODULE LEVEL CONVENIENCE FUNCTIONS/ User can import
#                                       without using the class

# ==========================================================================================================


def fetch_sequences(
    version,
    taxon_ids=None,
    proteome_id=None,
    go_id=None,
    pfam_id=None,
    db_config=None,
):
    """
    One-call helper: connect → query → disconnect → return records.


    Example
    -------
    >>> from get_reference_uniprot_set_lib import fetch_sequences
    >>> records = fetch_sequences("2026_01", taxon_ids=[9606, 10090])
    >>> print(len(records), "proteins retrieved")
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins(
            version=version,
            taxon_ids=taxon_ids,
            proteome_id=proteome_id,
            go_id=go_id,
            pfam_id=pfam_id,
        )


def fetch_sequences_by_hmm_hit(
    version,
    hmm_query,
    evalue_cutoff=None,
    taxon_ids=None,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → disconnect → return records.

    Retrieves sequences for all proteins with a hit to the given Pfam HMM profile.
    Results are deduplicated — one sequence per protein regardless of domain copy count.
    Suitable as a direct input to alignment and tree-building pipelines.


    Examples
    --------
    >>> from get_reference_uniprot_set_lib import fetch_sequences_by_hmm_hit
    >>> # All Homeodomain proteins across all taxa
    >>> records = fetch_sequences_by_hmm_hit("2026_01", "Homeodomain")
    >>> print(len(records), "proteins retrieved")

    >>> # Human + mouse kinases, strict E-value, pipe into BioPython
    >>> from get_reference_uniprot_set_lib import UniProtRetriever, get_db_config
    >>> records = fetch_sequences_by_hmm_hit(
    ...     "2026_01", "PF00069", evalue_cutoff=1e-10, taxon_ids=[9606, 10090]
    ... )
    >>> with UniProtRetriever(get_db_config()) as db:
    ...     seqrecords = db.to_biopython(records)
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins_by_hmm_hit(
            version=version,
            hmm_query=hmm_query,
            evalue_cutoff=evalue_cutoff,
            taxon_ids=taxon_ids,
        )


def fetch_fasta_string(
    version,
    taxon_ids=None,
    proteome_id=None,
    go_id=None,
    pfam_id=None,
    db_config=None,
):
    """
    One-call helper: connect → query → return a FASTA string.

    Useful for sending directly into Bio.SeqIO.parse() without writing a file.

    Example
    -------
    >>> import io
    >>> from Bio import SeqIO
    >>> from get_reference_uniprot_set_lib import fetch_fasta_string
    >>> fasta = fetch_fasta_string("2026_01", taxon_ids=[9606])
    >>> seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        records = db.get_proteins(
            version=version,
            taxon_ids=taxon_ids,
            proteome_id=proteome_id,
            go_id=go_id,
            pfam_id=pfam_id,
        )
        return db.to_fasta_string(records)


def fetch_fasta_string_by_hmm_hit(
    version,
    hmm_query,
    evalue_cutoff=None,
    taxon_ids=None,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → return a FASTA string.

    Useful for piping HMM-filtered sequences directly into Bio.SeqIO.parse()
    without writing a file.


    Example
    -------
    >>> import io
    >>> from Bio import SeqIO
    >>> from get_reference_uniprot_set_lib import fetch_fasta_string_by_hmm_hit
    >>> fasta = fetch_fasta_string_by_hmm_hit("2026_01", "Homeodomain", taxon_ids=[9606])
    >>> seqs = list(SeqIO.parse(io.StringIO(fasta), "fasta"))
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        records = db.get_proteins_by_hmm_hit(
            version=version,
            hmm_query=hmm_query,
            evalue_cutoff=evalue_cutoff,
            taxon_ids=taxon_ids,
        )
        return db.to_fasta_string(records)


def fetch_domains_by_accession(
    version,
    accessions,
    evalue_cutoff=None,
    chunk_size=5000,
    db_config=None,
):
    """
    One-call helper: connect → query hmm_search_results → disconnect →
    return all domain hits for the given accession(s).

    Returns one dict per domain occurrence (not one per protein). Use
    this when you need domain coordinates, scores, and HMM accessions
    rather than sequences.


    Examples
    --------
    >>> from get_reference_uniprot_set_lib import fetch_domains_by_accession
    >>> domains = fetch_domains_by_accession("2026_01", "P04637")
    >>> for d in domains:
    ...     print(d["hmm_name"], d["hmm_accession"], d["ali_from"], d["ali_to"])

    >>> # Works with any size list — chunked automatically
    >>> domains = fetch_domains_by_accession("2026_01", my_list_of_1000_accessions)
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_domains_by_accession(
            version=version,
            accessions=accessions,
            evalue_cutoff=evalue_cutoff,
            chunk_size=chunk_size,
        )


def fetch_domains_by_go(version, go_id, evalue_cutoff=None, db_config=None):
    """Return all HMM profiles found in proteins annotated with a GO term."""
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_domains_by_go(
            version=version, go_id=go_id, evalue_cutoff=evalue_cutoff
        )


def fetch_sequences_by_accession(
    version,
    accessions,
    chunk_size=5000,
    db_config=None,
):
    """
    One-call helper: connect → query proteins + sequences → disconnect →
    return full protein records for the given accession(s).


    Examples
    --------
    >>> from get_reference_uniprot_set_lib import fetch_sequences_by_accession
    >>> records = fetch_sequences_by_accession("2026_01", "P04637")
    >>> for r in records:
    ...     print(r["accession"], r["organism"])

    >>> # Works with any size list — chunked automatically
    >>> records = fetch_sequences_by_accession("2026_01", my_list_of_accessions)
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_proteins_by_accession(
            version=version,
            accessions=accessions,
            chunk_size=chunk_size,
        )


def fetch_presence_absence_matrix(
    version,
    pfam_queries,
    taxon_ids=None,
    evalue_cutoff=None,
    db_config=None,
):
    """
    One-call helper: connect → build presence/absence matrix → disconnect.


    Returns
    -------
    list[dict]  keys: taxon_id, scientific_name, hmm_name, hmm_accession,
                      hmm_type, protein_count
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_presence_absence_matrix(
            version=version,
            pfam_queries=pfam_queries,
            taxon_ids=taxon_ids,
            evalue_cutoff=evalue_cutoff,
        )


def fetch_accessions_for_cell(
    version,
    pfam_query,
    taxon_id,
    evalue_cutoff=None,
    db_config=None,
):
    """
    One-call helper: connect → get accessions for one matrix cell → disconnect.


    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_accessions_for_cell(
            version=version,
            pfam_query=pfam_query,
            taxon_id=taxon_id,
            evalue_cutoff=evalue_cutoff,
        )


def fetch_subprofile_hits(
    version,
    accessions,
    evalue_cutoff=None,
    exclude_queries=None,
    db_config=None,
):
    """
    One-call helper for Step 2, Path A (sub-profile enrichment).

    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_subprofile_hits(
            version=version,
            accessions=accessions,
            evalue_cutoff=evalue_cutoff,
            exclude_queries=exclude_queries,
        )


def fetch_domain_architectures(
    version,
    accessions,
    evalue_cutoff=None,
    collapse_repeats=True,
    db_config=None,
):
    """
    One-call helper for Step 2, Path B (domain architecture patterns).



    Returns
    -------
    list[dict]  keys: architecture, arch_accessions, protein_count,
                      coverage, example_accessions
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_domain_architectures(
            version=version,
            accessions=accessions,
            evalue_cutoff=evalue_cutoff,
            collapse_repeats=collapse_repeats,
        )


def fetch_accessions_with_taxids(version, accessions, db_config=None):
    """
    One-call helper: connect --> bulk accession --> taxid lookup --> disconnect.
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_accessions_with_taxids(version, accessions)


def fetch_highres_profile(
    version,
    pfam_subclade_map,
    taxon_ids=None,
    binary=False,
    db_config=None,
):
    """
    One-call helper: connect --> build high-res phylogenetic profile --> disconnect.
    """
    config = db_config or get_db_config()
    with UniProtRetriever(config) as db:
        return db.get_highres_profile(
            version=version,
            pfam_subclade_map=pfam_subclade_map,
            taxon_ids=taxon_ids,
            binary=binary,
        )


# ==========================================================================================================

#                                                CLI

# ==========================================================================================================


def _build_parser():
    parser = argparse.ArgumentParser(
        description="UniProt Reference Set Retrieval Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "-version", required=True, help="UniProt version (e.g., 2026_01)"
    )
    parser.add_argument(
        "-taxonomy", nargs="+", type=int, help="One or more Taxonomy IDs"
    )
    parser.add_argument(
        "--proteome-id", help="Filter by Proteome ID (e.g., UP000005640)"
    )
    parser.add_argument("--go-id", help="Filter by GO ID (e.g., GO:0005634)")
    parser.add_argument("--pfam-id", help="Filter by Pfam ID (e.g., PF00870)")
    parser.add_argument(
        "--list-versions", action="store_true", help="List all available versions"
    )
    parser.add_argument(
        "--list-proteomes",
        action="store_true",
        help="List proteome IDs for this version",
    )
    parser.add_argument(
        "--output-dir", default=None, help="Directory for the output FASTA file"
    )
    parser.add_argument(
        "--hmm-name",
        default=None,
        help="Pfam HMM name or accession, e.g. Homeodomain or PF00046",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=None,
        help="Max full-sequence E-value (default: no filter for --domains, 1e-5 for --hmm-name)",
    )

    parser.add_argument(
        "--accession",
        nargs="+",
        metavar="ACCESSION",
        help="Fetch sequences for one or more protein accessions (e.g. P04637 P10275)",
    )

    parser.add_argument(
        "--domains",
        nargs="+",
        metavar="ACCESSION",
        help="Get all Pfam domains for one or more protein accessions (e.g. P04637 P10275)",
    )

    parser.add_argument(
        "--go-domains",
        metavar="GO_ID",
        help="Get all HMM Profiles found in proteins annotated with a GO term",
    )

    # Direct credentials (alternative to .env)
    parser.add_argument(
        "--host", default=None, help="Database host (overrides .env / default)"
    )
    parser.add_argument(
        "--user", default=None, help="Database user (overrides .env / default)"
    )
    parser.add_argument(
        "--password", default=None, help="Database password (overrides .env / default)"
    )
    parser.add_argument(
        "--database", default=None, help="Database name (overrides .env / default)"
    )

    return parser


def main():
    args = _build_parser().parse_args()

    if args.user and args.password is None:
        args.password = getpass.getpass(
            prompt=f"Password for {args.user}@{args.host or 'localhost'}: "
        )

    retriever = UniProtRetriever(
        get_db_config(
            host=args.host,
            user=args.user,
            password=args.password,
            database=args.database,
        )
    )

    try:
        retriever.connect()

        if args.list_versions:
            retriever.list_available_versions()
            return

        if args.list_proteomes:
            proteomes = retriever.get_proteome_ids(args.version)
            print(f"\nProteomes in version {args.version}:")
            for p_id in proteomes:
                print(f"  - {p_id}")
            return

        if args.accession:
            print(f"\n{'='*60}")
            print(f"Accession lookup — version: {args.version}")
            print(f"Accessions : {', '.join(args.accession)}")
            print(f"{'='*60}\n")
            records = retriever.get_proteins_by_accession(
                version=args.version,
                accessions=args.accession,
            )
            if not records:
                print(" No sequences found for the given accession(s).")
                return
            identifier = f"acc_{'_'.join(args.accession[:3])}"
            retriever.export_fasta(
                records,
                args.version,
                identifier,
                output_dir=args.output_dir,
            )
            print(f"\nSuccessfully retrieved {len(records):,} sequences.")
            return

        if args.domains:
            print(f"\n{'='*60}")
            print(f"Domain lookup — version: {args.version}")
            print(f"Accessions : {', '.join(args.domains)}")
            if args.evalue:
                print(f"E-value    : ≤ {args.evalue}")
            print(f"{'='*60}\n")
            domains = retriever.get_domains_by_accession(
                version=args.version,
                accessions=args.domains,
                evalue_cutoff=args.evalue,
            )
            if not domains:
                print(" No domains found for the given accession(s).")
                return
            print(
                f"{'Accession':<12} {'Domain':<25} {'Accession':<12} {'Start':>6} {'End':>6} {'E-value':>10} {'Score':>8}"
            )
            print("-" * 80)
            for d in domains:
                print(
                    f"{d['accession']:<12} {d['hmm_name']:<25} {d['hmm_accession']:<12} "
                    f"{d['ali_from']:>6} {d['ali_to']:>6} {d['full_evalue']:>10.2e} {d['full_score']:>8.1f}"
                )
            print(f"\nTotal domain hits: {len(domains)}")
            return

        if args.go_domains:
            print(f"\n{'='*60}")
            print(f"GO → Domain Profile lookup")
            print(f"GO Term  : {args.go_domains}")
            print(f"Version  : {args.version}")
            if args.evalue:
                print(f"E-value  : ≤ {args.evalue}")
            print(f"{'='*60}\n")
            results = retriever.get_domains_by_go(
                version=args.version,
                go_id=args.go_domains,
                evalue_cutoff=args.evalue,
            )
            if not results:
                print("No domain profiles found for this GO term.")
                return
            print(f"{'HMM Name':<25} {'Accession':<12} {'Protein Count':>15}")
            print("-" * 55)
            for r in results:
                print(
                    f"{r['hmm_name']:<25} {r['hmm_accession']:<12} {r['protein_count']:>15,}"
                )
            print(f"\nTotal profiles: {len(results)}")
            return

        print(f"\n{'='*60}")
        print(f"UniProt Reference Set Retrieval")
        print(f"Version : {args.version}")
        if args.taxonomy:
            print(f"Taxonomy: {args.taxonomy}")
        if args.proteome_id:
            print(f"Proteome: {args.proteome_id}")
        if args.go_id:
            print(f"GO Term : {args.go_id}")
        if args.pfam_id:
            print(f"Pfam ID : {args.pfam_id}")

        if args.hmm_name:
            print(f"HMM    : {args.hmm_name} (E-value ≤ {args.evalue})")
            print(f"{'='*60}\n")
            records = retriever.get_proteins_by_hmm_hit(
                version=args.version,
                hmm_query=args.hmm_name,
                evalue_cutoff=args.evalue,
                taxon_ids=args.taxonomy,
            )
            identifier = f"hmm_{args.hmm_name}"
        else:
            print(f"{'='*60}\n")
            records = retriever.get_proteins(
                version=args.version,
                taxon_ids=args.taxonomy,
                proteome_id=args.proteome_id,
                go_id=args.go_id,
                pfam_id=args.pfam_id,
            )
            identifier = "filtered_set"
            if args.proteome_id:
                identifier = args.proteome_id
            elif args.taxonomy:
                identifier = f"tax_{args.taxonomy[0]}"

        if not records:
            print("\n No matching data found.")
            return

        retriever.export_fasta(
            records,
            args.version,
            identifier,
            output_dir=args.output_dir,
        )
        print(f"\nSuccessfully retrieved {len(records):,} sequences.")

    except mysql.connector.Error as err:
        print(f"\n Database error: {err}")
        sys.exit(1)
    except Exception as e:
        print(f"\n Error: {e}")
        sys.exit(1)
    finally:
        retriever.close()


if __name__ == "__main__":
    main()
