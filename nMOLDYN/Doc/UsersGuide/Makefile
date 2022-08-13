pdf:
      mkdir PDF
	latex nMOLDYN_for_PDF.tex    
	makeindex nMOLDYN_for_PDF.acr -s nMOLDYN_for_PDF.ist -o nMOLDYN_for_PDF.acn
	latex nMOLDYN_for_PDF.tex
	dvipdf nMOLDYN_for_PDF.dvi PDF/nMOLDYN_ug.pdf

html:
      mkdir HTML
	latex nMOLDYN_for_HTML.tex
	makeindex nMOLDYN_for_HTML.acr -s nMOLDYN_for_HTML.ist -o nMOLDYN_for_HTML.acn
	latex nMOLDYN_for_HTML.tex
	cp Figures/logo.gif HTML/
	latex2html -split 8 -local_icons -dir HTML nMOLDYN_for_HTML.tex
	latex2html -split 8 -ps_images -local_icons -dir HTML nMOLDYN_for_HTML.tex
	awk '{gsub("<A NAME=\"tex2html_wrap16368\"", "", $$0); print > FILENAME}' HTML/index.html
	awk '{gsub(".*IMAGE gif.*", "<IMG WIDTH=500 HEIGHT=400 SRC=\"logo.gif\">", $$0); print > FILENAME}' HTML/index.html
