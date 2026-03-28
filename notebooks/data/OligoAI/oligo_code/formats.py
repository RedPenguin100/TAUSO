from pydantic import BaseModel
from pydantic import BaseModel, Field
from typing import List, Literal

class SQL(BaseModel):
    sql_command: str

class ColumnDescription(BaseModel):
    column_name: str
    description: str

class Table(BaseModel):
    table_description: str
    column_descriptions: list[ColumnDescription]
    csv: str

class Script(BaseModel):
    pyscript: str

class Modification(BaseModel):
    modification: str = Field(
        ...,
        title="Modification Name",
        description="Name of the modification. For sugar: MOE, cEt, LNA, OMe. For backbone: PS, PO"
    )
    type: Literal["sugar", "backbone"] = Field(
        ...,
        title="Modification Type",
        description="Type of the modification - either sugar or backbone"
    )
    positions: List[int] = Field(
        ...,
        title="Modification Positions",
        description="List of 1-based positions where this modification occurs. For backbone modifications, position N refers to the linkage between nucleotides N and N+1"
    )

class Chemistry(BaseModel):
    length: int = Field(
        ...,
        title="Sequence Length",
        description="Total length of the oligonucleotide sequence",
    )
    modifications: List[Modification] = Field(
        ...,
        title="Modifications",
        description="List of all modifications applied to the sequence. Any positions not specified use default chemistry (D for sugar, PO for backbone)"
    )

class RowChemistry(BaseModel):
    process_row_py: str

class Metadata(BaseModel):
    target_mrna: str
    target_gene: str
    cell_line: str
    cell_line_species: Literal["human", "mouse", "other"]
    dosage: float
    cells_per_well: int
    transfection_method: Literal["Electroporation", "Gymnosis", "Lipofection", "Other"]