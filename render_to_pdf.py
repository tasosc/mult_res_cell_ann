"""
Render report to PDF
"""
from io import BytesIO
from datetime import datetime
from os import PathLike
from fpdf import FPDF
import fpdf
from matplotlib.figure import Figure
from pandas import DataFrame

class RenderToPdf:
    """
    Render report to PDF
    """
    def __init__(self, vebosity_limit: int):
        self.pdf = FPDF()
        self.add_page()
        self.verbosity_limit = vebosity_limit
        self.pdf.set_font('Times', size=10)
        self.output_filename = "report.pdf"
    def set_frontpage(self, filename: str):
        text = f"""
        <section>
<h1>Multi cell annotation report</h1>
<p>
<b>Date</b>{datetime.now().strftime("%a %d %B %Y %H:%M:%S %Z")}
</p>
<p>
<b>sc RNA seq dataset</b>{filename}
</p>
</section>
"""
        self.pdf.write_html(text)
    def add_page(self):
        self.pdf.add_page()
    def start_title(self, title: str):
        self.pdf.add_page()
        self.pdf.write_html(f'<h2>{title}</h2>')
    
    def should_write(self, expected_verbosity: int):
        if (expected_verbosity > self.verbosity_limit):
            return False
        return True

    def write_text(self, text: str, expected_verbosity : int = 1):
        if not self.should_write(expected_verbosity):
            return
        self.pdf.write_html(f"<p>{text}</p>")
    def write_table(self, df: DataFrame, expected_verbosity : int = 2):
        if not self.should_write(expected_verbosity):
            return
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

    def write_fig(self, fig: Figure, title: str = None, expected_verbosity: int = 1):
        if not self.should_write(expected_verbosity):
            return
        if title:
            self.write_text(title, expected_verbosity)
        img_buf = BytesIO()
        fig.savefig(img_buf, format="png", dpi=200)
        self.pdf.image(img_buf, x= fpdf.enums.Align.C, y=40, w=self.pdf.epw - 40, keep_aspect_ratio=True)
        img_buf.close()

    def set_filename(self, filename: PathLike):
        self.output_filename=filename

    def close(self):
        self.pdf.output(self.output_filename)
