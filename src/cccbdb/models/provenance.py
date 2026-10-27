"""Data provenance and bibliographic reference models for scientific reproducibility."""

from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from typing import Optional

from pydantic import BaseModel as _PydanticBaseModel, Field
try:
    from pydantic import ConfigDict
except ImportError:
    def ConfigDict(**kwargs): return kwargs

class BaseModel(_PydanticBaseModel):
    if not hasattr(_PydanticBaseModel, "model_dump"):
        def model_dump(self, **kwargs):
            return self.dict(**kwargs)
    if not hasattr(_PydanticBaseModel, "model_dump_json"):
        def model_dump_json(self, **kwargs):
            return self.json(**kwargs)
    if not hasattr(_PydanticBaseModel, "model_copy"):
        def model_copy(self, **kwargs):
            return self.copy(**kwargs)


from cccbdb.config import PARSER_VERSION


class ReferenceCitation(BaseModel):
    """Literature citation or reference squib recorded in NIST CCCBDB."""

    model_config = ConfigDict(frozen=True)

    squib: str = Field(..., description="NIST CCCBDB reference identifier (e.g., '1979HUB/HER', 'CODATA')")
    citation: Optional[str] = Field(default=None, description="Full bibliographic reference string")
    doi: Optional[str] = Field(default=None, description="Digital Object Identifier (DOI) if available")


class DataProvenance(BaseModel):
    """Complete provenance record linking a parsed datum back to NIST CCCBDB."""

    model_config = ConfigDict(frozen=True)

    source_url: str = Field(..., description="Exact NIST CCCBDB URL from which the datum was extracted")
    retrieved_at: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),
        description="ISO-8601 UTC timestamp of retrieval",
    )
    parser_version: str = Field(default=PARSER_VERSION, description="Version of cccbdb-py parser")
    raw_content_sha256: Optional[str] = Field(
        default=None, description="SHA-256 hex digest of the raw HTML source page"
    )
    original_unit: Optional[str] = Field(
        default=None, description="Original physical unit reported by NIST CCCBDB"
    )
    reference_squib: Optional[str] = Field(
        default=None, description="Literature reference squib if experimental"
    )
    nist_release: str = Field(
        default="Release 22 (May 2022)", description="NIST CCCBDB database release identifier"
    )

    @classmethod
    def from_html(
        cls,
        source_url: str,
        html: str,
        original_unit: Optional[str] = None,
        reference_squib: Optional[str] = None,
    ) -> "DataProvenance":
        """Create a provenance record computing the SHA-256 digest of the source HTML."""
        digest = hashlib.sha256(html.encode("utf-8", errors="replace")).hexdigest()
        return cls(
            source_url=source_url,
            raw_content_sha256=digest,
            original_unit=original_unit,
            reference_squib=reference_squib,
        )
