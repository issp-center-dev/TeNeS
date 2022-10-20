ROOT_DIR=`pwd`
if [ -z "$(grep 'project(TeNeS' $ROOT_DIR/CMakeLists.txt 2>/dev/null)" ]; then
  echo "ERROR: current directory is not the root directory of the TeNeS codes"
  exit 1
fi

if [ ! -d $ROOT_DIR/.git ]; then
  echo "ERROR: this is not a git repository"
  exit 1
fi

res=0
type git-archive-all >/dev/null 2>&1 || res=1
if [ $res -eq 1 ]; then
  echo "ERROR: git-archive-all is not installed"
  exit 1
fi

cd $ROOT_DIR
rm -rf build-doc
mkdir build-doc
cd build-doc
cmake -DDocument=ON ../
for lang in ja en; do
  make doc-${lang}-pdf
  cp docs/sphinx/${lang}/pdf/TeNeS_${lang}.pdf ../
done
cd $ROOT_DIR

version=1.4-dev

git submodule update -i -r
git-archive-all \
  --prefix=TeNeS-${version} \
  --extra=TeNeS_ja.pdf \
  --extra=TeNeS_en.pdf \
  TeNeS-${version}.tar.gz
