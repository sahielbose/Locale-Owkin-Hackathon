"""Lane A analysis engine.

Every engine function takes the canonical AnnData and returns schema objects.
Do NOT put analysis logic anywhere else; the MCP layer only wraps this engine.

Submodules are imported directly (not re-exported here) so that Lane B can keep
importing each engine function lazily: a broken edit to one module isolates to the
tool that uses it instead of breaking every engine import.
"""
