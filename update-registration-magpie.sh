#!/bin/bash
# run this script with one or more arguments that are the directories that hold the code files for
# the objects you want to update.

dirs=$@

appname='MagpieApp'

objs='
AtomicDensityAux
MyTRIMDensityAux
MyTRIMElementEnergyAux
MyTRIMElementResultAux
SPPARKSAux
CoupledDefectAnnihilation
DefectAnnihilation
MyTRIMElementHeatSource
MyTRIMElementSource
MyTRIMDiracSource
MyTRIMDiracRun
MyTRIMElementRun
MyTRIMPKAInConeInfo
MyTRIMPKAInfo
MyTRIMRasterizer
NeutronicsSpectrumSamplerFission
NeutronicsSpectrumSamplerSN
PKAConstant
PKAFissionFragmentEmpirical
PKAFissionFragmentNeutronics
PKAFixedPointGenerator
PKAFunction
PKAGeneratorAlphaDecay
PKAGeneratorRecoil
PKAGun
SPPARKSUserObject
ElasticRecoil
InelasticRecoil
IsotopeRecoilRate
MultiAppNeutronicsSpectrumTransfer
MyTRIMMesh
MyTRIMDiracResult
MyTRIMPKAEnergyHistogram
MyTRIMPKAStatistics
PKAList
IsotopeRecoilRateSampler
...
'

actions='
'

for obj in $objs; do
    find $dirs | grep "/${obj}\.C$" | xargs perl -0777 -i'' -pe 's/(template ?<>\nInputParameters\n? ?validParams<(\w*)>\(\))/registerMooseObject\("'"$appname"'", \2\);\n\n\1/igs'
done

while read action task; do
    find $dirs | grep "/${action}\.C$" | xargs perl -0777 -i'' -pe 's/(template ?<>\nInputParameters\n? ?validParams<(\w*)>\(\))/registerMooseAction\("'"$appname"'", \2, "'"$task"'"\);\n\n\1/igs'
done <<< "$actions"
