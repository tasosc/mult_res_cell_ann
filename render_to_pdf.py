"""
Render report to PDF
"""
from io import BytesIO
from datetime import datetime
from os import PathLike
from fpdf import FPDF
from matplotlib.figure import Figure
from pandas import DataFrame

class RenderToPdf:
    """
    Render report to PDF
    """
    def __init__(self):
        self.pdf = FPDF()
        self.add_page()
        self.pdf.set_font('Times', size=10)
        self.output_filename = "report.pdf"
    def set_frontpage(self, filename: str):
        text = f"""
**Multi cell annotation report**

{datetime.now().strftime("%a %d %B %Y %H:%M:%S %Z")}

sc RNA seq dataset: {filename}
"""
        self.pdf.multi_cell(text=text,w=0.0, markdown=True)
    def add_page(self):
        self.pdf.add_page()
    def start_title(self, title: str):
        self.pdf.add_page()
        self.pdf.cell(text=f"**{title}**", markdown=True)
    def write_text(self, text: str):
        self.pdf.multi_cell(text=text, w=self.pdf.w, markdown=True)
    def write_table(self, df: DataFrame):
        df_str = df.applymap(str)
        columns = df.columns
        rows = df_str.values.tolist()

        with self.pdf.table(
            borders_layout="MINIMAL",
            cell_fill_color=200,  # grey
            cell_fill_mode="ROWS",
            line_height=self.pdf.font_size * 2,
            text_align="CENTER",
            width=180) as table:
            hader_row = table.row()
            for column in columns:
                hader_row.cell(column)
            for data_row in rows:
                row = table.row()
                for value in data_row:
                    row.cell(value)
    
    
    def write_fig(self, fig: Figure, title: str = None, alt_text: str = None):
       img_buf = BytesIO()
       fig.savefig(img_buf, format="svg")
       self.pdf.image(img_buf, title=title, alt_text=alt_text)
       img_buf.close()

    def set_filename(self, filename: PathLike):
        self.output_filename=filename
    
    def close(self):
        self.pdf.output(self.output_filename)
    








