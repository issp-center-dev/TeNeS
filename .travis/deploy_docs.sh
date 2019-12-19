# This is a pull request, finish.
if [ "_$TRAVIS_PULL_REQUEST" != "_false" ] ;then
  echo "This is a pull request, do nothing."
  exit 0;
fi

if [ "_$TRAVIS_BRANCH" == "_master" ]; then
  echo "This is the master branch, deploy docs."
elif [ -n "$TRAVIS_TAG" ]; then
  echo "This is a versioned tag, deploy docs."
else
  echo "Do nothing."
  exit 0
fi

openssl aes-256-cbc -K $encrypted_d3b6bca1dfe3_key -iv $encrypted_d3b6bca1dfe3_iv -in ${ROOTDIR}/.travis/id_rsa.enc -out ~/.ssh/id_rsa -d
chmod 600 ~/.ssh/id_rsa
echo -e "Host github.com\n\tStrictHostKeyChecking no\n" >> ~/.ssh/config

cd ${ROOTDIR}
git clone git@github.com:${TRAVIS_REPO_SLUG} tenes-repo
cd tenes-repo
git checkout gh-pages
cd manual
if [ "_${TRAVIS_BRANCH}" == "_master" ]; then
  rm -rf master
  mkdir -p master/ja
  mkdir -p master/en
  cp -r ${ROOTDIR}/build/docs/sphinx/ja/html master/ja
  cp -r ${ROOTDIR}/build/docs/sphinx/en/html master/en
  git add master
elif [ -n ${TRAVIS_TAG} ]; then
  rm -rf ${TRAVIS_TAG}
  mkdir -p ${TRAVIS_TAG}/ja
  mkdir -p ${TRAVIS_TAG}/en
  cp -r ${ROOTDIR}/build/docs/sphinx/ja/html ${TRAVIS_TAG}/ja
  cp -r ${ROOTDIR}/build/docs/sphinx/en/html ${TRAVIS_TAG}/en
  git add ${TRAVIS_TAG}
else
  echo "The deploy script failed to solve where to install documents. The script has some mistake."
  echo "\$TRAVIS_BRANCH: $TRAVIS_BRANCH"
  echo "\$TRAVIS_TAG: $TRAVIS_TAG"
  echo "\$TRAVIS_PULL_REQUEST: $TRAVIS_PULL_REQUEST"
  exit 1
fi

git config --global user.name "TeNeS developers"
git config --global user.email "tenes-dev@issp.u-tokyo.ac.jp"
git commit -m "Update by TravisCI (#${TRAVIS_BUILD_NUMBER})"
ST=$?
if [ $ST == 0 ]; then
  git push origin gh-pages:gh-pages --follow-tags > /dev/null 2>&1
fi

