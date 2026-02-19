# CCCBDB Server Status Report

*Generated: 2026-02-18 22:12:22*
*Source: 8 fetch run(s)*

## Overview

This report logs all HTTP requests made during data archival from the CCCBDB
website (https://cccbdb.nist.gov), a NIST-maintained ASP.NET application behind
Cloudflare CDN. Transient errors (timeouts, 503s, most 500s) are expected during
bulk fetching and resolve with retries — the error rate dropped from ~11% on the
initial daytime pass to 0% on subsequent off-peak retries.

## Aggregate Statistics

| Metric | Value |
|--------|-------|
| Total HTTP requests | 26,322 |
| Successful | 25,162 (95.6%) |
| Failed | 1,160 (4.4%) |
| Avg response time (success) | 489 ms |

## Fetch Runs

| Run Start | Run End | Requests | Errors | Error Rate |
|-----------|---------|----------|--------|------------|
| 2026-02-16T09:54:51 | 2026-02-16T15:30:54 | 5,100 | 561 | 11.0% |
| 2026-02-16T17:27:27 | 2026-02-16T17:27:35 | 6 | 0 | 0.0% |
| 2026-02-16T17:28:31 | 2026-02-16T17:28:39 | 6 | 0 | 0.0% |
| 2026-02-16T17:46:16 | 2026-02-16T19:23:00 | 3,961 | 62 | 1.6% |
| 2026-02-16T19:30:31 | 2026-02-16T21:46:59 | 5,039 | 3 | 0.1% |
| 2026-02-18T21:28:44 | 2026-02-18T21:28:52 | 20 | 0 | 0.0% |
| 2026-02-18T21:29:25 | 2026-02-18T21:38:31 | 1,260 | 0 | 0.0% |
|  |  | 10,930 | 534 | 4.9% |

## Error Types

| Error Type | Count | % of Errors | Description |
|------------|-------|-------------|-------------|
| 500 | 965 | 83.2% | Internal Server Error — transient under load; permanent for some geom3x species |
| timeout | 107 | 9.2% | Read timeout (>30s) — transient, more common during peak hours |
| 503 | 87 | 7.5% | Service Unavailable — transient Cloudflare/server load, resolves on retry |
| Remote end closed connection without response | 1 | 0.1% |  |

## Error Rate by Page Type

| Page Type | Total Requests | Successes | Errors | Error Rate | Top Error |
|-----------|---------------|-----------|--------|------------|-----------|
| expgeom2x | 9,963 | 9,873 | 90 | 0.9% | timeout |
| geom2x | 6,558 | 6,351 | 207 | 3.2% | 500 |
| geom3x | 3,482 | 2,673 | 809 | 23.2% | 500 |
| spin2x | 4,061 | 4,034 | 27 | 0.7% | 503 |
| energy2x | 2,258 | 2,231 | 27 | 1.2% | 500 |

## Error Timeline (by hour)

| Hour | Successes | Errors | Error Rate |
|------|-----------|--------|------------|
| 2026-02-16T09 | 60 | 8 | 11.8% |
| 2026-02-16T10 | 798 | 96 | 10.7% |
| 2026-02-16T11 | 939 | 78 | 7.7% |
| 2026-02-16T12 | 831 | 102 | 10.9% |
| 2026-02-16T13 | 753 | 110 | 12.7% |
| 2026-02-16T14 | 535 | 127 | 19.2% |
| 2026-02-16T15 | 623 | 40 | 6.0% |
| 2026-02-16T17 | 617 | 0 | 0.0% |
| 2026-02-16T18 | 2362 | 14 | 0.6% |
| 2026-02-16T19 | 2032 | 48 | 2.3% |
| 2026-02-16T20 | 2261 | 0 | 0.0% |
| 2026-02-16T21 | 1675 | 3 | 0.2% |
| 2026-02-18T21 | 1280 | 0 | 0.0% |

## Fetch Behavior Notes

The error timeline above shows a clear pattern: the initial daytime fetch had
~11% errors, but successive retries during off-peak hours (evenings/weekends)
drove the rate to 0%. All errors except a small number of permanent `geom3x`
failures were transient and caused by server load, not by bugs in the server
or the client.

### Transient errors resolve with retries

HTTP 500, 503, and timeout errors are load-dependent. The `--missing-only` flag
retries only species still missing data, and running during off-peak US Eastern
hours yields near-zero error rates. Session failures (stale server-side state)
are automatically detected and recovered by `--heal`.

### Permanent `geom3x` 500 errors

A subset of species consistently return HTTP 500 on `geom3x.asp` (calculated
geometry detail pages). The server cannot compute certain method/basis
combinations for these molecules. These are the only genuinely permanent
failures — experimental geometry, spin, and energy data for the same species
typically succeed.

### Technical notes for developers

- **Session requirement**: CCCBDB uses server-side ASP.NET sessions. The initial
  `expgeom2x.asp?casno=X&charge=Y` request establishes the molecule context;
  subsequent pages (`spin2x.asp`, `energy2x.asp`, etc.) rely on the session.
  The mirror always fetches `expgeom2x.asp` first.
- **Use `urllib`, not `curl`**: Cloudflare's challenge-platform token (set via
  JavaScript) must be propagated. Python's `HTTPCookieProcessor` handles this;
  `curl` with exported cookies does not.
- **HTML tag inconsistency**: `expgeom2x.asp` uses unclosed `<TD>` tags while
  `geom3x.asp` uses closed `</TD>` tags. Parsers use `(?:</TD>)?` to handle both.

## Most Error-Prone Species (Top 20)

| Species | Error Count |
|---------|-------------|
| C2H3N | 7 |
| CH3NO | 6 |
| C3H6O | 5 |
| C3H8S | 5 |
| C2H5N | 4 |
| C2H5O | 4 |
| CH2N2 | 4 |
| CH3N2+ | 4 |
| FNO2 | 4 |
| FNS | 4 |
| C2HF(+1)_2713099 | 4 |
| CHFO_1493023 | 4 |
| CHFO(+1)_1493023 | 4 |
| BF2H_13709836 | 4 |
| CHF2_2670135 | 4 |
| CHF2(+1)_2670135 | 4 |
| F2HN_10405273 | 4 |
| CH2F2_75105 | 4 |
| C2H3F_75025 | 4 |
| C2H3F(+1)_75025 | 4 |

## Re-fetch Strategy

1. **Retry with `--missing-only`** — fills gaps left by transient errors.
2. **Run during off-peak hours** (evenings/weekends US Eastern) — error rates
   drop from ~11% to near 0%.
3. **Run `--heal`** to detect and recover session failures (stale server-side
   state saved as bad HTML). Scans, clears, re-fetches, loops until convergence.
4. **Accept permanent `geom3x` 500s** — these species genuinely lack computed
   geometry on the server. All other data types succeed for them.
5. **2-second rate limit** between requests is hardcoded. Do not reduce it.

---
*Report covers 8 fetch run(s) with 26,322 total requests.*
