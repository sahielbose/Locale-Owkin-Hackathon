"""Locale: an MCP server that reasons about where cells sit in tumor tissue.

Subpackages:
  schema      shared Pydantic data contract (the coordination artifact)
  data        build the canonical AnnData from extracted CSVs (Lane A)
  engine      spatial analysis (graph, enrichment, niches, outcome, validate) (Lane A)
  mcp_server  MCP tools wrapping the engine (Lane B)
  viz         tissue-map payload + widget (Lane C)
"""

__version__ = "0.1.0"
