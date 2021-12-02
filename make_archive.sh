rm -rf build-doc
mkdir build-doc
cd build-doc
cmake -DDocument=ON ../
for lang in ja en; do
  make doc-${lang}-pdf
  cp docs/sphinx/${lang}/pdf/TeNeS_${lang}.pdf ../
done
cd ../

version=1.2.0

git submodule update -i -r
git-archive-all \
  --prefix=TeNeS-${version} \
  --extra=TeNeS_ja.pdf \
  --extra=TeNeS_en.pdf \
  TeNeS-${version}.tar.gz
