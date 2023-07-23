# build html version of documentation
# after buiding, you can open the index.html in build/html/
sphinx-build -a -b html . build/html

# build pdf version of documentation
#sphinx-build -a -b latex . build/latex
#make -C build/latex all-pdf
#cp build/latex/WannierTools.pdf .
