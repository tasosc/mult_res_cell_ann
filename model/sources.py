from pydantic import BaseModel


class Sources(BaseModel):
    names: list[str]
