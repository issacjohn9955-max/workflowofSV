"""Render the Chinese manual from its maintained UTF-8 text source."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from PIL import ImageFont

ROOT=Path(__file__).resolve().parents[1]
FONT=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')


def lines_for(text, measure, limit):
    lines=[]
    line=''
    for char in text:
        if line and measure.getlength(line+char)>limit:
            lines.append(line)
            line=char.lstrip()
        else:
            line+=char
    return lines+[line]


def render(font,out):
    face=FontProperties(fname=str(font))
    measure=ImageFont.truetype(str(font),100)
    pages=(ROOT/'docs/manual_zh.txt').read_text().strip().split('\n---PAGE---\n')
    with plt.rc_context({'pdf.fonttype':3}), PdfPages(out,metadata={'Title':'SV PC3、PC5、PC10 工作流手册'}) as pdf:
        for index,chapter in enumerate(pages,1):
            title,*body=chapter.splitlines()
            fig=plt.figure(figsize=(8.27,11.69),facecolor='white')
            fig.text(.085,.945,'SV–ENVIRONMENT ASSOCIATION  /  METHODS AND OPERATIONS',
                     color='#386C8C',fontsize=9,va='top')
            fig.text(.085,.897,title,color='#273746',fontsize=20,fontproperties=face,va='top')
            fig.add_artist(plt.Line2D([.085,.915],[.852,.852],transform=fig.transFigure,color='#A9B2BB',lw=.7))
            y=.821
            for paragraph in body:
                if paragraph.startswith('§ '):
                    fig.text(.085,y,paragraph[2:],color='#386C8C',fontsize=11.4,fontproperties=face,va='top')
                    y-=.036
                elif paragraph.startswith('https://'):
                    fig.text(.085,y,paragraph,color='#386C8C',fontsize=9,va='top',url=paragraph)
                    y-=.035
                else:
                    lines=lines_for(paragraph,measure,4280)
                    for line in lines:
                        fig.text(.085,y,line,color='#273746',fontsize=11.15,fontproperties=face,va='top')
                        y-=.027
                    y-=.012
            if y<.075:
                raise ValueError(f'Manual page {index} exceeds available space.')
            fig.add_artist(plt.Line2D([.085,.915],[.068,.068],transform=fig.transFigure,color='#A9B2BB',lw=.55))
            fig.text(.085,.044,'LFMM2 · Partial RDA · PC3 / PC5 / PC10',color='#386C8C',fontsize=8)
            fig.text(.915,.044,f'{index:02d} / {len(pages):02d}',ha='right',color='#386C8C',fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--font',type=Path,default=FONT)
    parser.add_argument('--out',type=Path,default=ROOT/'docs/SV_PC_LFMM2_RDA_操作与原理.pdf')
    args=parser.parse_args()
    render(args.font,args.out)
