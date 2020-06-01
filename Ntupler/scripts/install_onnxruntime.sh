#!/bin/bash
set -e
set -x

if [[ -z "$CMSSW_BASE" ]]; then
	echo '$CMSSW_BASE not defined -- run cmsenv first!'
fi

curdir=`pwd`
INSTALL_DIR="$CMSSW_BASE/ORT_INSTALL"

[[ -d $INSTALL_DIR ]] && rm -r $INSTALL_DIR
mkdir $INSTALL_DIR && cd $INSTALL_DIR

curl -s http://hqu.web.cern.ch/hqu/tools/onnxruntime/onnxruntime-1.2.0_ppc_58f53966d.tar.gz -o ort.tgz
tar -xf ort.tgz && rm ort.tgz

sed -i 's@_ONNXRUNTIME_INSTALL_PATH_@'"$INSTALL_DIR"'@g' $INSTALL_DIR/onnxruntime.xml
scram setup $INSTALL_DIR/onnxruntime.xml

cd $CMSSW_BASE/src
git cms-addpkg PhysicsTools/ONNXRuntime
cd $CMSSW_BASE/src/PhysicsTools/ONNXRuntime
scram b clean && scram b -j16

cp -a --remove-destination $INSTALL_DIR/lib/libonnxruntime.so* $CMSSW_BASE/external/$SCRAM_ARCH/lib
rm -r $INSTALL_DIR/lib

cd $curdir
