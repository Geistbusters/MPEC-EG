root=`pwd`
cd ..;chmod og+rx $root;cd $root
### for my mac m1 laptop
#f1="mpif90"
#f2="mpif90 -fallow-argument-mismatch"
#f3=$root"/lib/lib"
#opt=-O2
#opt=-Og
#optno=-Og
#uselib=1
#cat <<@ > ~/.prempec
#alias runmpec='mpirun -n 6 $root/mpec'
#@
### for NAS facility Pleaides
f1="ifort -axAVX -mcmodel=medium -xSSE4.1"
f2=$f1
f3="-lmpi -mkl -shared-intel"
#opt="-O0 -g -traceback -check all -ftrapuv"
opt=-O2
optno=-O0
uselib=0
cat <<@ > ~/.prempec
eval \`/usr/bin/modulecmd bash load comp-intel/2020.4.304 mpi-hpe/mpt\`
PATH=\$PATH:/u/scicon/tools/bin
alias runmpec='mpiexec mbind.x $root/mpec'
@
chmod og+rx ~/.prempec
. ~/.prempec
### ddi1
ddi=1
### ddi4
#ddi=4
### library dgemm
#dgemm=0
### debug dgemm
dgemm=1
### end of choices
if [ -e $root/test ]
then
echo "test exits!"
else
echo "test needs to be untared"
tar -xf test.tar
cd test
chmod og+r *
cd ..
chmod og+rx test
fi
s=$root/source
if [ -e $s ]
then
echo "source exits!"
else
echo "source needs to be untared"
tar -xf source.tar
fi
/bin/rm mpec
if [ $ddi -eq 1 ]
then
sed s/cddi1// < $s/ddi.template > $s/ddi.f
sed s/cddi1// < $s/common.mympi.template > $s/common.mympi
else
sed s/cddi4// < $s/ddi.template > $s/ddi.f
sed s/cddi4// < $s/common.mympi.template > $s/common.mympi
fi
if [ $dgemm -eq 0 ]
then
cp $s/start.template $s/start.f
else
sed s/cdgm// < $s/start.template > $s/start.f
fi
if [ -e $root/obj ]
then
echo "obj exists! "
else
echo "obj needs to be created"
mkdir $root/obj
fi
if [ -e $root/data ]
then
echo "data exists! "
else
echo "data needs to be untared"
tar -xf $root/data.tar
cd $root/data
chmod og+r *
cd ..
chmod og+xr data
fi
if [ $uselib -ne 0 ]
then
if [ -e $root/lib ]
then
echo "lib exits!"
else
echo "creating lib"
tar xf $root/lib.tar
fi
cd $root/lib
cat > makefilelib<<!
lib: dasum.o daxpy.o dcopy.o ddot.o dgemm.o dgemv.o dger.o dlacpy.o dlae2.o dlaebz.o dlaev2.o dlagtf.o dlagts.o dlamch.o dlansy.o dlapy2.o dlarf.o dlarfb.o dlarfg.o \
dlarft.o dlarnv.o dlartg.o dlaruv.o dlasr.o dlassq.o dlatrd.o dlazro.o dnrm2.o dorg2l.o dorg2r.o dorgql.o dorgqr.o dorgtr.o dorm2l.o dorm2r.o dormql.o dormqr.o \
dormtr.o drot.o drotg.o dscal.o dstebz.o dstein.o dsteqr.o dsterf.o dswap.o dsyevx.o dsymv.o dsyr2.o dsyr2k.o dsytd2.o dsytrd.o dtrmm.o dtrmv.o idamax.o ilaenv.o \
lsame.o xerbla.o dstevx.o dlanst.o dgesvd.o dbdsqr.o dgebrd.o dgelqf.o dgeqrf.o dlange.o dlascl.o dlaset.o dorgbr.o dorglq.o dormbr.o dgebd2.o dgelq2.o dgeqr2.o \
disnan.o dlabrd.o dlas2.o dlasq1.o dlasv2.o dorgl2.o dormlq.o dlaisnan.o dlasq2.o dlasrt.o dorml2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o
	ar -q lib dasum.o daxpy.o dcopy.o ddot.o dgemm.o dgemv.o dger.o dlacpy.o dlae2.o dlaebz.o dlaev2.o dlagtf.o dlagts.o dlamch.o dlansy.o dlapy2.o dlarf.o dlarfb.o dlarfg.o \
dlarft.o dlarnv.o dlartg.o dlaruv.o dlasr.o dlassq.o dlatrd.o dlazro.o dnrm2.o dorg2l.o dorg2r.o dorgql.o dorgqr.o dorgtr.o dorm2l.o dorm2r.o dormql.o dormqr.o \
dormtr.o drot.o drotg.o dscal.o dstebz.o dstein.o dsteqr.o dsterf.o dswap.o dsyevx.o dsymv.o dsyr2.o dsyr2k.o dsytd2.o dsytrd.o dtrmm.o dtrmv.o idamax.o ilaenv.o \
lsame.o xerbla.o dstevx.o dlanst.o dgesvd.o dbdsqr.o dgebrd.o dgelqf.o dgeqrf.o dlange.o dlascl.o dlaset.o dorgbr.o dorglq.o dormbr.o dgebd2.o dgelq2.o dgeqr2.o \
disnan.o dlabrd.o dlas2.o dlasq1.o dlasv2.o dorgl2.o dormlq.o dlaisnan.o dlasq2.o dlasrt.o dorml2.o dlasq3.o dlasq4.o dlasq5.o dlasq6.o
dlasq4.o : dlasq4.f ; $f2 -c $opt dlasq4.f
dlasq5.o : dlasq5.f ; $f2 -c $opt dlasq5.f
dlasq6.o : dlasq6.f ; $f2 -c $opt dlasq6.f
dlasq3.o : dlasq3.f ; $f2 -c $opt dlasq3.f
dlaisnan.o : dlaisnan.f ; $f2 -c $opt dlaisnan.f
dlasq2.o : dlasq2.f ; $f2 -c $opt dlasq2.f
dlasrt.o : dlasrt.f ; $f2 -c $opt dlasrt.f
dorml2.o : dorml2.f ; $f2 -c $opt dorml2.f
dgebd2.o : dgebd2.f ; $f2 -c $opt dgebd2.f
dgelq2.o : dgelq2.f ; $f2 -c $opt dgelq2.f
dgeqr2.o : dgeqr2.f ; $f2 -c $opt dgeqr2.f
disnan.o : disnan.f ; $f2 -c $opt disnan.f
dlabrd.o : dlabrd.f ; $f2 -c $opt dlabrd.f
dlas2.o : dlas2.f ; $f2 -c $opt dlas2.f
dlasq1.o : dlasq1.f ; $f2 -c $opt dlasq1.f
dlasv2.o : dlasv2.f ; $f2 -c $opt dlasv2.f
dorgl2.o : dorgl2.f ; $f2 -c $opt dorgl2.f
dormlq.o : dormlq.f ; $f2 -c $opt dormlq.f
dbdsqr.o : dbdsqr.f ; $f2 -c $opt dbdsqr.f
dgebrd.o : dgebrd.f ; $f2 -c $opt dgebrd.f
dgelqf.o : dgelqf.f ; $f2 -c $opt dgelqf.f
dgeqrf.o : dgeqrf.f ; $f2 -c $opt dgeqrf.f
dlange.o : dlange.f ; $f2 -c $opt dlange.f
dlascl.o : dlascl.f ; $f2 -c $opt dlascl.f
dlaset.o : dlaset.f ; $f2 -c $opt dlaset.f
dorgbr.o : dorgbr.f ; $f2 -c $opt dorgbr.f
dorglq.o : dorglq.f ; $f2 -c $opt dorglq.f
dormbr.o : dormbr.f ; $f2 -c $opt dormbr.f
dgesvd.o : dgesvd.f ; $f2 -c $opt dgesvd.f
dstevx.o : dstevx.f ; $f2 -c $opt dstevx.f
dlanst.o : dlanst.f ; $f2 -c $opt dlanst.f
dasum.o : dasum.f ; $f2 -c $opt dasum.f
daxpy.o : daxpy.f ; $f2 -c $opt daxpy.f
dcopy.o : dcopy.f ; $f2 -c $opt dcopy.f
ddot.o : ddot.f ; $f2 -c $opt ddot.f
dgemm.o : dgemm.f ; $f2 -c $opt dgemm.f
dgemv.o : dgemv.f ; $f2 -c $opt dgemv.f
dger.o : dger.f ; $f2 -c $opt dger.f
dlacpy.o : dlacpy.f ; $f2 -c $opt dlacpy.f
dlae2.o : dlae2.f ; $f2 -c $opt dlae2.f
dlaebz.o : dlaebz.f ; $f2 -c $opt dlaebz.f
dlaev2.o : dlaev2.f ; $f2 -c $opt dlaev2.f
dlagtf.o : dlagtf.f ; $f2 -c $opt dlagtf.f
dlagts.o : dlagts.f ; $f2 -c $opt dlagts.f
dlamch.o : dlamch.f ; $f2 -c $opt dlamch.f
dlansy.o : dlansy.f ; $f2 -c $opt dlansy.f
dlapy2.o : dlapy2.f ; $f2 -c $opt dlapy2.f
dlarf.o : dlarf.f ; $f2 -c $opt dlarf.f
dlarfb.o : dlarfb.f ; $f2 -c $opt dlarfb.f
dlarfg.o : dlarfg.f ; $f2 -c $opt dlarfg.f
dlarft.o : dlarft.f ; $f2 -c $opt dlarft.f
dlarnv.o : dlarnv.f ; $f2 -c $opt dlarnv.f
dlartg.o : dlartg.f ; $f2 -c $opt dlartg.f
dlaruv.o : dlaruv.f ; $f2 -c $opt dlaruv.f
dlasr.o : dlasr.f ; $f2 -c $opt dlasr.f
dlassq.o : dlassq.f ; $f2 -c $opt dlassq.f
dlatrd.o : dlatrd.f ; $f2 -c $opt dlatrd.f
dlazro.o : dlazro.f ; $f2 -c $opt dlazro.f
dnrm2.o : dnrm2.f ; $f2 -c $opt dnrm2.f
dorg2l.o : dorg2l.f ; $f2 -c $opt dorg2l.f
dorg2r.o : dorg2r.f ; $f2 -c $opt dorg2r.f
dorgql.o : dorgql.f ; $f2 -c $opt dorgql.f
dorgqr.o : dorgqr.f ; $f2 -c $opt dorgqr.f
dorgtr.o : dorgtr.f ; $f2 -c $opt dorgtr.f
dorm2l.o : dorm2l.f ; $f2 -c $opt dorm2l.f
dorm2r.o : dorm2r.f ; $f2 -c $opt dorm2r.f
dormql.o : dormql.f ; $f2 -c $opt dormql.f
dormqr.o : dormqr.f ; $f2 -c $opt dormqr.f
dormtr.o : dormtr.f ; $f2 -c $opt dormtr.f
drot.o : drot.f ; $f2 -c $opt drot.f
drotg.o : drotg.f ; $f2 -c $opt drotg.f
dscal.o : dscal.f ; $f2 -c $opt dscal.f
dstebz.o : dstebz.f ; $f2 -c $opt dstebz.f
dstein.o : dstein.f ; $f2 -c $opt dstein.f
dsteqr.o : dsteqr.f ; $f2 -c $opt dsteqr.f
dsterf.o : dsterf.f ; $f2 -c $opt dsterf.f
dswap.o : dswap.f ; $f2 -c $opt dswap.f
dsyevx.o : dsyevx.f ; $f2 -c $opt dsyevx.f
dsymv.o : dsymv.f ; $f2 -c $opt dsymv.f
dsyr2.o : dsyr2.f ; $f2 -c $opt dsyr2.f
dsyr2k.o : dsyr2k.f ; $f2 -c $opt dsyr2k.f
dsytd2.o : dsytd2.f ; $f2 -c $opt dsytd2.f
dsytrd.o : dsytrd.f ; $f2 -c $opt dsytrd.f
dtrmm.o : dtrmm.f ; $f2 -c $opt dtrmm.f
dtrmv.o : dtrmv.f ; $f2 -c $opt dtrmv.f
idamax.o : idamax.f ; $f2 -c $opt idamax.f
ilaenv.o : ilaenv.f ; $f2 -c $opt ilaenv.f
lsame.o : lsame.f ; $f2 -c $opt lsame.f
xerbla.o : xerbla.f ; $f2 -c $opt xerbla.f
!
make -f makefilelib lib
fi
cd $root/obj
cat > makefile<<!
mpec: start.o mkmem.o addcomma4.o basisz.o delim.o ddi.o cartsfromv.o copyto.o runit.o enough.o delimb.o njsymb16.o genr.o getm.o getxtr.o gopt.o \
intin.o loadr2.o loadr.o mostpop.o parseterm.o prntm2.o lsqfit2.o randd8.o spher.o gaussqnetlib.o lusolvd.o plmv.o parahf.o \
buildcasgrad.o buildhesscas.o cannon.o cas0.o contractg.o dhf.o cas1.o cas1b.o cas1bnona1.o cgvec.o cgvec2.o chj.o werup.o updateg.o srtsymc.o \
srtsym.o square.o spindig.o sortoutsym.o sortdet.o sgen.o savetonormal.o savetops.o secoi.o setbqn.o pshamdbl.o pshamdnona1.o pshampp.o pshamppnona1.o \
psl.o psspin.o pshamd.o prtocc.o propcas.o prop.o printat.o sogen.o solin.o sospher2.o sospher3.o stepwisecsf.o wrot.o prtsigs.o putwf.o relgammas.o \
wrot1.o wrot2.o wrot2c.o precsf.o printa.o wfetch.o parap.o phouse.o phsdet.o parah0.o parajkfromh.o paramp2.o parap42.o phasv.o paraeri42cf.o \
paraeri4c.o paraeri.o parah042c.o parah04c.o orbsymspace.o paracis.o restin.o sldvrc.o ssdvrc.o opto.o oplist.o taylor.o onei.o onep.o \
moldenf.o mpprnt2.o mrci.o msl.o nattrac.o reloadvs.o reordergv.o savewf.o tofro.o trans2e.o updatexc.o updatexuc.o vnormx.o mofromao.o prephss.o \
prephdd.o prtit.o puttoiprop.o recoverwf.o make1dm.o make4xdm.o makeguess.o maptoold.o mkdoub.o mkdoubuc.o mksing.o ucdoubop.o loaddoubb.o loaddoubc.o \
loadsinga.o loadsingb.o m12gen.o m1gen.o int1tobit.o intcsf.o intcsfcas.o intsum.o loaddouba.o pslzzcsfcas.o pythag.o r2sb.o strip2ps.o \
stuffintohiicsf.o updateiicsf.o updateq.o lzziicsfcas.o mov2ps.o prodn.o prtcsfvec.o prtcsfveccas.o pshamcsf.o pshamcsfcas.o pslzzcsf.o stuffintohii2.o \
ilimts.o int1tobit2.o int1tobit3.o lzziicsf.o uncompxu.o hiicsf.o hiicsfcas.o hsingpart.o hsscsf.o ifind2.o hcis.o hcsi.o hcss.o hddcsf.o \
hdoubpart.o hdoubpartuc.o hessycmo.o hessyccas.o hsscsf12.o hcdduc.o hcdi.o hcdiuc.o hcds.o hcdsuc.o hcid.o hciduc.o hddcsf12b.o xtimesn.o hccsfcas.o \
hccsfd.o hccsfrilr.o hcdd4v.o hcdddiag.o hcddjk.o hcdduc1.o hcdduc2.o xtimesn2.o hc1cd.o hc1cdbl.o hc1cnona1.o hc1cnona1d.o hccsf.o mxmt.o paraerik.o \
sandtd.o vdmo2so.o getbas0.o getwf.o gotoqs.o graborb.o grado.o gradcas.o gtmom.o gtmomso.o hc1c.o parajkfromhd2c.o psioppsi.o rordr.o sym4o.o tofrob.o \
transder.o transder1.o trianglez.o parajkfromhd1j.o parajkfromhd2b.o testme.o unmakegd.o unmakegs.o parah0d2.o parah0grad.o parajkfromhd0.o parajkfromhd1.o \
sotest2.o sotest2c.o sotest.o onedints.o orbder.o orbdercas.o paraerid2b.o paraerid2c.o paraerid.o paraeridj.o int4copy.o makegd.o makegs.o mtimesx.o nattracd.o \
oned2ints.o oneider.o restid.o restind.o genwf.o hcbk.o hccsfbk.o hccsfd1.o hcddjkd1.o hcisd1.o hcssd1.o hcdsbk.o hcidbk.o hcidbkuc.o hcisbk.o hcsdbk.o hcsibk.o \
hcssbk.o wnorm.o get4int.o hcddjkbk.o hcdducbk.o hcdibk.o hcdibkuc.o hcdsd1.o wnorm1.o wnorm2.o wnorm2c.o genryd.o genmatn.o genmatn2.o singnm1.o genmat.o \
genhsoa2.o genfcnp.o psisopsi.o genfcn2.o genf0.o genergy.o gendertype.o gencsf2b.o gencsf3.o hcidbk4.o hcisbk4.o hcsdbk4.o hcsibk4.o hcssbk4.o i4to8copy.o \
paraerikbk4.o testsome.o vdmo2sobk4.o getcup.o hccsfbk4.o hcddjkbk4.o hcdibk4.o hcdsbk4.o spinloop.o spinloop1.o gen12nona1.o gencsf1.o gencsf2.o gencup.o \
uncompxut.o gcode.o gcsfps.o gcsfpslz.o gen12.o epdvec.o erir.o erird.o fiddleh.o fillrup.o gandc.o gandc4.o gandcr.o dynwtr.o \
dumpbas0.o dumpbas2.o dumpbas.o dumpdoubc.o dumpdoubuc.o dumph.o dumpsing.o dosingnona1.o dosingnona1d.o dotdvecs.o dotvall.o dovr.o dumb.o dumpbas1.o \
doab2dbl.o dodoub.o dodoubd.o dodoubdbl.o dodoubnona1.o dodoubnona1d.o dosing.o dosingd.o dosingdbl.o doab11.o doab1.o doab1d.o doab1dbl.o doab1non.o \
doab1nond.o doab2.o doab2d.o doab21.o diagx.o diagy.o distit8.o distit.o do4v.o derofxor.o dgerid.o diagdy.o diaghii.o diaghiid.o diaghiidbl.o diaglzz.o \
derh01b.o derh0.o derid.o derofxor1.o diaghii1.o diaglzz1.o deconvert2torc.o deconvert12torc.o denmx.o denmx12.o der3part.o deramat.o derh01.o hccsfd12.o \
hcsid12.o corenota1.o corelz2.o countcc.o countem2.o countfcn.o cupo1rp.o cupo21.o dcannon.o dcbit.o cont2ecsf.o cont2ecsfa.o cont2ecsfb.o convert12tor.o \
convert12torc.o convert1tor.o convert2tor.o convert2torc.o core.o consolidate1.o matchc.o orthovrb.o compxtimes.o consolidate3.o consolidate2.o matchcl.o \
check1.o check2.o check3.o checkc.o checkps.o cmpvcsf.o cnt12.o cnt12nona1.o compxtimes2.o bodcpart.o buildder3s.o buildhess.o buildhessdx.o c2eprop.o \
cal1int.o cart2spher.o cfileget.o atomll.o avgr.o mapv2v.o pdump.o addcomma8.o addqtop.o bdens.o hccsfdx.o genmatn3.o prehcddjk.o cmpvcsf2.o cleb2.o \
getint.o getintcas.o hcdsuc1.o hcdsucbk.o hcdsuc1ds.o hcsducbk.o hcdsuc1sd.o igetint.o paraeri4ct.o assqn.o dsortdws.o timescomp.o move2min.o minimize.o \
ntrans.o hcssd12.o eina.o run.o compcomp.o prntm2r.o mtimesh.o mtimesh2.o getgss.o tdendet.o mtimesi.o dorb2e.o dorbmix.o dorbd4o.o dorbd1x.o dorbdj.o \
dorbdk.o genivder.o foldr.o intcsfder.o diagxd.o precder.o hcdi12.o hcds12.o hcds12c.o hcddjkd12.o den14oc.o dorbd3x.o chc4v.o dorbd4x.o trans3xden.o dorb4v.o \
abtrans3.o tt4v.o paraeridd.o der4v.o postden.o
	 $f1 -o mpec start.o addcomma4.o basisz.o delim.o ddi.o mkmem.o cartsfromv.o copyto.o runit.o enough.o delimb.o njsymb16.o genr.o getm.o getxtr.o gopt.o \
intin.o loadr2.o loadr.o mostpop.o parseterm.o prntm2.o lsqfit2.o randd8.o spher.o gaussqnetlib.o lusolvd.o plmv.o parahf.o \
buildcasgrad.o buildhesscas.o cannon.o cas0.o contractg.o dhf.o cas1.o cas1b.o cas1bnona1.o cgvec.o cgvec2.o chj.o werup.o updateg.o srtsymc.o \
srtsym.o square.o spindig.o sortoutsym.o sortdet.o sgen.o savetonormal.o savetops.o secoi.o setbqn.o pshamdbl.o pshamdnona1.o pshampp.o pshamppnona1.o \
psl.o psspin.o pshamd.o prtocc.o propcas.o prop.o printat.o sogen.o solin.o sospher2.o sospher3.o stepwisecsf.o wrot.o prtsigs.o putwf.o relgammas.o \
wrot1.o wrot2.o wrot2c.o precsf.o printa.o wfetch.o parap.o phouse.o phsdet.o parah0.o parajkfromh.o paramp2.o parap42.o phasv.o paraeri42cf.o \
paraeri4c.o paraeri.o parah042c.o parah04c.o orbsymspace.o paracis.o restin.o sldvrc.o ssdvrc.o opto.o oplist.o taylor.o onei.o onep.o \
moldenf.o mpprnt2.o mrci.o msl.o nattrac.o reloadvs.o reordergv.o savewf.o tofro.o trans2e.o updatexc.o updatexuc.o vnormx.o mofromao.o prephss.o \
prephdd.o prtit.o puttoiprop.o recoverwf.o make1dm.o make4xdm.o makeguess.o maptoold.o mkdoub.o mkdoubuc.o mksing.o ucdoubop.o loaddoubb.o loaddoubc.o \
loadsinga.o loadsingb.o m12gen.o m1gen.o int1tobit.o intcsf.o intcsfcas.o intsum.o loaddouba.o pslzzcsfcas.o pythag.o r2sb.o strip2ps.o \
stuffintohiicsf.o updateiicsf.o updateq.o lzziicsfcas.o mov2ps.o prodn.o prtcsfvec.o prtcsfveccas.o pshamcsf.o pshamcsfcas.o pslzzcsf.o stuffintohii2.o \
ilimts.o int1tobit2.o int1tobit3.o lzziicsf.o uncompxu.o hiicsf.o hiicsfcas.o hsingpart.o hsscsf.o ifind2.o hcis.o hcsi.o hcss.o hddcsf.o \
hdoubpart.o hdoubpartuc.o hessycmo.o hessyccas.o hsscsf12.o hcdduc.o hcdi.o hcdiuc.o hcds.o hcdsuc.o hcid.o hciduc.o hddcsf12b.o xtimesn.o hccsfcas.o \
hccsfd.o hccsfrilr.o hcdd4v.o hcdddiag.o hcddjk.o hcdduc1.o hcdduc2.o xtimesn2.o hc1cd.o hc1cdbl.o hc1cnona1.o hc1cnona1d.o hccsf.o mxmt.o paraerik.o \
sandtd.o vdmo2so.o getbas0.o getwf.o gotoqs.o graborb.o grado.o gradcas.o gtmom.o gtmomso.o hc1c.o parajkfromhd2c.o psioppsi.o rordr.o sym4o.o tofrob.o \
transder.o transder1.o trianglez.o parajkfromhd1j.o parajkfromhd2b.o testme.o unmakegd.o unmakegs.o parah0d2.o parah0grad.o parajkfromhd0.o parajkfromhd1.o \
sotest2.o sotest2c.o sotest.o onedints.o orbder.o orbdercas.o paraerid2b.o paraerid2c.o paraerid.o paraeridj.o int4copy.o makegd.o makegs.o mtimesx.o nattracd.o \
oned2ints.o oneider.o restid.o restind.o genwf.o hcbk.o hccsfbk.o hccsfd1.o hcddjkd1.o hcisd1.o hcssd1.o hcdsbk.o hcidbk.o hcidbkuc.o hcisbk.o hcsdbk.o hcsibk.o \
hcssbk.o wnorm.o get4int.o hcddjkbk.o hcdducbk.o hcdibk.o hcdibkuc.o hcdsd1.o wnorm1.o wnorm2.o wnorm2c.o genryd.o genmatn.o genmatn2.o singnm1.o genmat.o \
genhsoa2.o genfcnp.o psisopsi.o genfcn2.o genf0.o genergy.o gendertype.o gencsf2b.o gencsf3.o hcidbk4.o hcisbk4.o hcsdbk4.o hcsibk4.o hcssbk4.o i4to8copy.o \
paraerikbk4.o testsome.o vdmo2sobk4.o getcup.o hccsfbk4.o hcddjkbk4.o hcdibk4.o hcdsbk4.o spinloop.o spinloop1.o gen12nona1.o gencsf1.o gencsf2.o gencup.o \
uncompxut.o gcode.o gcsfps.o gcsfpslz.o gen12.o epdvec.o erir.o erird.o fiddleh.o fillrup.o gandc.o gandc4.o gandcr.o dynwtr.o \
dumpbas0.o dumpbas2.o dumpbas.o dumpdoubc.o dumpdoubuc.o dumph.o dumpsing.o dosingnona1.o dosingnona1d.o dotdvecs.o dotvall.o dovr.o dumb.o dumpbas1.o \
doab2dbl.o dodoub.o dodoubd.o dodoubdbl.o dodoubnona1.o dodoubnona1d.o dosing.o dosingd.o dosingdbl.o doab11.o doab1.o doab1d.o doab1dbl.o doab1non.o \
doab1nond.o doab2.o doab2d.o doab21.o diagx.o diagy.o distit8.o distit.o do4v.o derofxor.o dgerid.o diagdy.o diaghii.o diaghiid.o diaghiidbl.o diaglzz.o \
derh01b.o derh0.o derid.o derofxor1.o diaghii1.o diaglzz1.o deconvert2torc.o deconvert12torc.o denmx.o denmx12.o der3part.o deramat.o derh01.o hccsfd12.o \
hcsid12.o corenota1.o corelz2.o countcc.o countem2.o countfcn.o cupo1rp.o cupo21.o dcannon.o dcbit.o cont2ecsf.o cont2ecsfa.o cont2ecsfb.o convert12tor.o \
convert12torc.o convert1tor.o convert2tor.o convert2torc.o core.o consolidate1.o matchc.o orthovrb.o compxtimes.o consolidate3.o consolidate2.o matchcl.o \
check1.o check2.o check3.o checkc.o checkps.o cmpvcsf.o cnt12.o cnt12nona1.o compxtimes2.o bodcpart.o buildder3s.o buildhess.o buildhessdx.o c2eprop.o \
cal1int.o cart2spher.o cfileget.o atomll.o avgr.o mapv2v.o pdump.o addcomma8.o addqtop.o bdens.o hccsfdx.o genmatn3.o prehcddjk.o cmpvcsf2.o cleb2.o \
getint.o getintcas.o hcdsuc1.o hcdsucbk.o hcdsuc1ds.o hcsducbk.o hcdsuc1sd.o igetint.o paraeri4ct.o assqn.o dsortdws.o timescomp.o move2min.o minimize.o \
ntrans.o hcssd12.o eina.o run.o compcomp.o prntm2r.o mtimesh.o mtimesh2.o getgss.o tdendet.o mtimesi.o dorb2e.o dorbmix.o dorbd4o.o dorbd1x.o dorbdj.o \
dorbdk.o genivder.o foldr.o intcsfder.o diagxd.o precder.o hcdi12.o hcds12.o hcds12c.o hcddjkd12.o den14oc.o dorbd3x.o dorbd4x.o trans3xden.o chc4v.o dorb4v.o \
abtrans3.o tt4v.o paraeridd.o der4v.o postden.o $f3
start.o : $s/start.f $s/common.store ; $f2 -c $opt $s/start.f
run.o : $s/run.f $s/common.basis $s/common.input $s/common.mrci $s/common.store; $f2 -c $opt $s/run.f
eina.o : $s/eina.f ; $f2 -c $opt $s/eina.f
abtrans3.o : $s/abtrans3.f $s/common.store ; $f2 -c $opt $s/abtrans3.f
tt4v.o : $s/tt4v.f $s/common.store ; $f2 -c $opt $s/tt4v.f
postden.o : $s/postden.f $s/common.store ; $f2 -c $opt $s/postden.f
precder.o : $s/precder.f $s/common.store ; $f2 -c $opt $s/precder.f
hcdi12.o : $s/hcdi12.f $s/common.store ; $f2 -c $opt $s/hcdi12.f
hcds12.o : $s/hcds12.f $s/common.store ; $f2 -c $opt $s/hcds12.f
hcddjkd12.o : $s/hcddjkd12.f $s/common.store ; $f2 -c $opt $s/hcddjkd12.f
den14oc.o : $s/den14oc.f $s/common.store ; $f2 -c $opt $s/den14oc.f
hcds12c.o : $s/hcds12c.f $s/common.store ; $f2 -c $opt $s/hcds12c.f
foldr.o : $s/foldr.f $s/common.store ; $f2 -c $opt $s/foldr.f
diagxd.o : $s/diagxd.f $s/common.store $s/common.hf; $f2 -c $opt $s/diagxd.f
intcsfder.o : $s/intcsfder.f $s/common.store $s/common.print ; $f2 -c $opt $s/intcsfder.f
dorb2e.o : $s/dorb2e.f $s/common.store ; $f2 -c $opt $s/dorb2e.f
dorbmix.o : $s/dorbmix.f $s/common.store ; $f2 -c $opt $s/dorbmix.f
dorbd4o.o : $s/dorbd4o.f $s/common.store ; $f2 -c $opt $s/dorbd4o.f
dorbd1x.o : $s/dorbd1x.f $s/common.store ; $f2 -c $opt $s/dorbd1x.f
dorbd3x.o : $s/dorbd3x.f $s/common.store ; $f2 -c $opt $s/dorbd3x.f
dorbd4x.o : $s/dorbd4x.f $s/common.store ; $f2 -c $opt $s/dorbd4x.f
trans3xden.o : $s/trans3xden.f $s/common.store ; $f2 -c $opt $s/trans3xden.f
dorb4v.o : $s/dorb4v.f $s/common.store ; $f2 -c $opt $s/dorb4v.f
chc4v.o : $s/chc4v.f $s/common.store ; $f2 -c $opt $s/chc4v.f
dorbdj.o : $s/dorbdj.f $s/common.store ; $f2 -c $opt $s/dorbdj.f
dorbdk.o : $s/dorbdk.f $s/common.store ; $f2 -c $opt $s/dorbdk.f
genivder.o : $s/genivder.f $s/common.basis $s/common.store ; $f2 -c $opt $s/genivder.f
mtimesi.o : $s/mtimesi.f $s/common.store ; $f2 -c $opt $s/mtimesi.f
mtimesh.o : $s/mtimesh.f ; $f2 -c $opt $s/mtimesh.f
mtimesh2.o : $s/mtimesh2.f $s/common.store ; $f2 -c $opt $s/mtimesh2.f
tdendet.o : $s/tdendet.f $s/common.store ; $f2 -c $opt $s/tdendet.f
getgss.o : $s/getgss.f $s/common.store ; $f2 -c $opt $s/getgss.f
prntm2r.o : $s/prntm2r.f ; $f2 -c $opt $s/prntm2r.f
compcomp.o : $s/compcomp.f ; $f2 -c $opt $s/compcomp.f
ntrans.o : $s/ntrans.f $s/common.store $s/common.basis $s/common.input ; $f2 -c $opt $s/ntrans.f                                                                        
hcssd12.o : $s/hcssd12.f $s/common.store ; $f2 -c $opt $s/hcssd12.f
move2min.o : $s/move2min.f $s/common.store $s/common.basis $s/common.input $s/common.mrci ; $f2 -c $opt $s/move2min.f
minimize.o : $s/minimize.f ; $f2 -c $opt $s/minimize.f
bdens.o : $s/bdens.f $s/common.store ; $f2 -c $opt $s/bdens.f
assqn.o : $s/assqn.f $s/common.store ; $f2 -c $opt $s/assqn.f
dsortdws.o : $s/dsortdws.f ; $f2 -c $opt $s/dsortdws.f
timescomp.o : $s/timescomp.f $s/common.store ; $f2 -c $opt $s/timescomp.f
getint.o : $s/getint.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/getint.f
igetint.o : $s/igetint.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/igetint.f
getintcas.o : $s/getintcas.f $s/common.store ; $f2 -c $opt $s/getintcas.f
hcdsuc1.o : $s/hcdsuc1.f $s/common.store ; $f2 -c $opt $s/hcdsuc1.f
hcdsucbk.o : $s/hcdsucbk.f $s/common.store ; $f2 -c $optno $s/hcdsucbk.f
hcdsuc1ds.o : $s/hcdsuc1ds.f $s/common.store ; $f2 -c $opt $s/hcdsuc1ds.f
hcsducbk.o : $s/hcsducbk.f $s/common.store ; $f2 -c $opt $s/hcsducbk.f
hcdsuc1sd.o : $s/hcdsuc1sd.f $s/common.store ; $f2 -c $opt $s/hcdsuc1sd.f
cleb2.o : $s/cleb2.f ; $f2 -c $opt $s/cleb2.f
hccsfdx.o : $s/hccsfdx.f $s/common.store $s/common.print ; $f2 -c $opt $s/hccsfdx.f
runit.o : $s/runit.f $s/common.print $s/common.mrci $s/common.cas $s/common.hf $s/common.store $s/common.rys $s/common.basis $s/common.spher $s/common.input $s/mpec.date ; $f2 -c $opt $s/runit.f
addcomma4.o : $s/addcomma4.f ; $f2 -c $opt $s/addcomma4.f
basisz.o : $s/basisz.f $s/cas.par $s/mrci.par $s/common.hf $s/common.cas $s/common.store $s/common.basis $s/common.input $s/common.spher $s/common.mrci $s/common.print ; $f2 -c $optno $s/basisz.f
delim.o : $s/delim.f ; $f2 -c $opt $s/delim.f
ddi.o : $s/ddi.f $s/common.store $s/common.mympi ; $f2 -c $opt $s/ddi.f
mkmem.o : $s/mkmem.f90 ; $f2 -c $s/mkmem.f90
cartsfromv.o : $s/cartsfromv.f $s/common.basis; $f2 -c $opt $s/cartsfromv.f
copyto.o : $s/copyto.f $s/common.store $s/common.basis $s/common.input $s/common.print ; $f2 -c $opt $s/copyto.f
enough.o : $s/enough.f ; $f2 -c $opt $s/enough.f
delimb.o : $s/delimb.f ; $f2 -c $opt $s/delimb.f
njsymb16.o : $s/njsymb16.f ; $f2 -c $opt $s/njsymb16.f
genr.o : $s/genr.f ; $f2 -c $opt $s/genr.f
getm.o : $s/getm.f ; $f2 -c $opt $s/getm.f
getxtr.o : $s/getxtr.f ; $f2 -c $opt $s/getxtr.f
gopt.o : $s/gopt.f $s/common.basis $s/common.input $s/common.store ; $f2 -c $opt $s/gopt.f
intin.o : $s/intin.f ; $f2 -c $opt $s/intin.f
loadr2.o : $s/loadr2.f $s/common.basis $s/common.input $s/common.store $s/common.print ; $f2 -c $opt $s/loadr2.f
loadr.o : $s/loadr.f $s/common.basis $s/common.input $s/common.store ; $f2 -c $opt $s/loadr.f
mostpop.o : $s/mostpop.f ; $f2 -c $opt $s/mostpop.f
parseterm.o : $s/parseterm.f ; $f2 -c $opt $s/parseterm.f
prntm2.o : $s/prntm2.f ; $f2 -c $opt $s/prntm2.f
lsqfit2.o : $s/lsqfit2.f ; $f2 -c $opt $s/lsqfit2.f
spher.o : $s/spher.f $s/common.spher $s/common.store ; $f2 -c $opt $s/spher.f
randd8.o : $s/randd8.f ; $f2 -c $opt $s/randd8.f
gaussqnetlib.o : $s/gaussqnetlib.f ; $f2 -c $opt $s/gaussqnetlib.f
lusolvd.o : $s/lusolvd.f ; $f2 -c $opt $s/lusolvd.f
plmv.o : $s/plmv.f ; $f2 -c $opt $s/plmv.f
parahf.o : $s/parahf.f $s/common.hf $s/common.cas $s/common.rys $s/common.spher $s/common.store $s/common.print $s/common.basis ; $f2 -c $opt $s/parahf.f
buildcasgrad.o : $s/buildcasgrad.f $s/common.store $s/common.hf ; $f2 -c $opt $s/buildcasgrad.f
buildhesscas.o : $s/buildhesscas.f $s/common.hf $s/common.store ; $f2 -c $opt $s/buildhesscas.f
cannon.o : $s/cannon.f $s/common.hf $s/common.store $s/common.print ; $f2 -c $opt $s/cannon.f
cas0.o : $s/cas0.f $s/common.store $s/common.hf $s/common.cas $s/common.print ; $f2 -c $opt $s/cas0.f
contractg.o : $s/contractg.f $s/common.store $s/common.hf $s/common.print ; $f2 -c $opt $s/contractg.f
dhf.o : $s/dhf.f $s/common.store ; $f2 -c $opt $s/dhf.f
cas1.o : $s/cas1.f $s/common.store $s/common.hf $s/common.cas ; $f2 -c $opt $s/cas1.f
cas1b.o : $s/cas1b.f $s/common.store $s/common.hf $s/common.cas ; $f2 -c $opt $s/cas1b.f
cas1bnona1.o : $s/cas1bnona1.f $s/common.store $s/common.hf $s/common.cas ; $f2 -c $opt $s/cas1bnona1.f
cgvec.o : $s/cgvec.f $s/common.store ; $f2 -c $opt $s/cgvec.f
cgvec2.o : $s/cgvec2.f $s/common.store ; $f2 -c $opt $s/cgvec2.f
chj.o : $s/chj.f $s/common.store $s/common.basis $s/common.input ; $f2 -c $opt $s/chj.f
werup.o : $s/werup.f $s/common.store $s/common.hf ; $f2 -c $opt $s/werup.f
updateg.o : $s/updateg.f $s/common.store $s/common.hf $s/common.print ; $f2 -c $opt $s/updateg.f
srtsymc.o : $s/srtsymc.f ; $f2 -c $opt $s/srtsymc.f
srtsym.o : $s/srtsym.f ; $f2 -c $opt $s/srtsym.f
square.o : $s/square.f ; $f2 -c $opt $s/square.f
spindig.o : $s/spindig.f $s/common.cas $s/common.store ; $f2 -c $opt $s/spindig.f
sortoutsym.o : $s/sortoutsym.f $s/common.store $s/common.print ; $f2 -c $opt $s/sortoutsym.f
sortdet.o : $s/sortdet.f $s/common.hf $s/common.cas $s/common.store ; $f2 -c $opt $s/sortdet.f
sgen.o : $s/sgen.f ; $f2 -c $opt $s/sgen.f
savetonormal.o : $s/savetonormal.f $s/common.store ; $f2 -c $opt $s/savetonormal.f
savetops.o : $s/savetops.f $s/common.store ; $f2 -c $opt $s/savetops.f
secoi.o : $s/secoi.f $s/common.store $s/common.hf ; $f2 -c $opt $s/secoi.f
setbqn.o : $s/setbqn.f $s/common.store $s/common.spher ; $f2 -c $opt $s/setbqn.f
pshamdbl.o : $s/pshamdbl.f $s/common.cas $s/common.store ; $f2 -c $opt $s/pshamdbl.f
pshamdnona1.o : $s/pshamdnona1.f $s/common.cas $s/common.store ; $f2 -c $opt $s/pshamdnona1.f
pshampp.o : $s/pshampp.f $s/common.cas $s/common.store ; $f2 -c $opt $s/pshampp.f
pshamppnona1.o : $s/pshamppnona1.f $s/common.cas $s/common.store ; $f2 -c $opt $s/pshamppnona1.f
psl.o : $s/psl.f $s/common.cas $s/common.store ; $f2 -c $opt $s/psl.f
psspin.o : $s/psspin.f $s/common.cas $s/common.store ; $f2 -c $opt $s/psspin.f
pshamd.o : $s/pshamd.f $s/common.cas $s/common.store ; $f2 -c $opt $s/pshamd.f
prtocc.o : $s/prtocc.f $s/common.store ; $f2 -c $opt $s/prtocc.f
propcas.o : $s/propcas.f $s/common.store $s/common.mrci $s/common.hf ; $f2 -c $opt $s/propcas.f
prop.o : $s/prop.f $s/common.store $s/common.mrci $s/common.hf $s/common.basis ; $f2 -c $opt $s/prop.f
printat.o : $s/printat.f ; $f2 -c $opt $s/printat.f
sogen.o : $s/sogen.f $s/common.store ; $f2 -c $opt $s/sogen.f
solin.o : $s/solin.f $s/common.store $s/common.print ; $f2 -c $opt $s/solin.f
sospher2.o : $s/sospher2.f $s/common.store $s/common.print ; $f2 -c $opt $s/sospher2.f
sospher3.o : $s/sospher3.f $s/common.store ; $f2 -c $opt $s/sospher3.f
stepwisecsf.o : $s/stepwisecsf.f $s/common.print $s/common.store ; $f2 -c $opt $s/stepwisecsf.f
wrot.o : $s/wrot.f $s/common.store ; $f2 -c $opt $s/wrot.f
prtsigs.o : $s/prtsigs.f ; $f2 -c $opt $s/prtsigs.f
putwf.o : $s/putwf.f $s/common.store $s/common.basis ; $f2 -c $opt $s/putwf.f
relgammas.o : $s/relgammas.f $s/common.store ; $f2 -c $opt $s/relgammas.f
wrot1.o : $s/wrot1.f $s/common.store ; $f2 -c $opt $s/wrot1.f
wrot2.o : $s/wrot2.f $s/common.store ; $f2 -c $opt $s/wrot2.f
wrot2c.o : $s/wrot2c.f $s/common.store ; $f2 -c $opt $s/wrot2c.f
precsf.o : $s/precsf.f $s/common.cas $s/common.store ; $f2 -c $opt $s/precsf.f
printa.o : $s/printa.f ; $f2 -c $opt $s/printa.f
wfetch.o : $s/wfetch.f $s/common.store ; $f2 -c $opt $s/wfetch.f
parap.o : $s/parap.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parap.f
phouse.o : $s/phouse.f $s/common.store ; $f2 -c $opt $s/phouse.f
phsdet.o : $s/phsdet.f ; $f2 -c $opt $s/phsdet.f
parah0.o : $s/parah0.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parah0.f
parajkfromh.o : $s/parajkfromh.f $s/common.store $s/common.hf $s/common.print ; $f2 -c $opt $s/parajkfromh.f
paramp2.o : $s/paramp2.f $s/common.print $s/common.store $s/common.hf ; $f2 -c $opt $s/paramp2.f
parap42.o : $s/parap42.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parap42.f
phasv.o : $s/phasv.f ; $f2 -c $opt $s/phasv.f
paraeri42cf.o : $s/paraeri42cf.f $s/common.hf $s/common.store $s/common.spher $s/common.print $s/common.basis ; $f2 -c $opt $s/paraeri42cf.f
paraeri4c.o : $s/paraeri4c.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/paraeri4c.f
paraeri4ct.o : $s/paraeri4ct.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/paraeri4ct.f
paraeri.o : $s/paraeri.f $s/common.hf $s/common.store $s/common.basis $s/common.print ; $f2 -c $opt $s/paraeri.f
der4v.o : $s/der4v.f $s/common.store ; $f2 -c $opt $s/der4v.f
paraeridd.o : $s/paraeridd.f $s/common.hf $s/common.store $s/common.basis $s/common.print ; $f2 -c $opt $s/paraeridd.f
parah042c.o : $s/parah042c.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parah042c.f
parah04c.o : $s/parah04c.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parah04c.f
orbsymspace.o : $s/orbsymspace.f $s/common.store ; $f2 -c $opt $s/orbsymspace.f
paracis.o : $s/paracis.f $s/common.store $s/common.hf ; $f2 -c $opt $s/paracis.f
restin.o : $s/restin.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/restin.f
sldvrc.o : $s/sldvrc.f $s/common.store ; $f2 -c $opt $s/sldvrc.f
ssdvrc.o : $s/ssdvrc.f $s/common.store ; $f2 -c $opt $s/ssdvrc.f
opto.o : $s/opto.f $s/common.store $s/common.hf $s/common.print ; $f2 -c $opt $s/opto.f
oplist.o : $s/oplist.f $s/common.basis $s/common.input ; $f2 -c $opt $s/oplist.f
taylor.o : $s/taylor.f $s/common.store ; $f2 -c $opt $s/taylor.f
onei.o : $s/onei.f $s/common.store $s/common.spher ; $f2 -c $opt $s/onei.f
onep.o : $s/onep.f $s/common.store $s/common.spher ; $f2 -c $optno $s/onep.f
moldenf.o : $s/moldenf.f $s/common.store ; $f2 -c $opt $s/moldenf.f
mpprnt2.o : $s/mpprnt2.f ; $f2 -c $opt $s/mpprnt2.f
mrci.o : $s/mrci.f $s/common.store $s/common.mrci $s/common.hf $s/common.print ; $f2 -c $opt $s/mrci.f
msl.o : $s/msl.f ; $f2 -c $opt $s/msl.f
nattrac.o : $s/nattrac.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/nattrac.f
reloadvs.o : $s/reloadvs.f $s/common.store ; $f2 -c $opt $s/reloadvs.f
reordergv.o : $s/reordergv.f $s/common.store ; $f2 -c $opt $s/reordergv.f
savewf.o : $s/savewf.f $s/common.store ; $f2 -c $opt $s/savewf.f
tofro.o : $s/tofro.f $s/common.store ; $f2 -c $opt $s/tofro.f
trans2e.o : $s/trans2e.f $s/common.store $s/common.hf ; $f2 -c $opt $s/trans2e.f
updatexc.o : $s/updatexc.f $s/common.basis $s/common.input $s/common.store $s/common.print ; $f2 -c $opt $s/updatexc.f
updatexuc.o : $s/updatexuc.f $s/common.basis $s/common.input $s/common.store $s/common.print ; $f2 -c $opt $s/updatexuc.f
vnormx.o : $s/vnormx.f ; $f2 -c $opt $s/vnormx.f
mofromao.o : $s/mofromao.f $s/common.hf $s/common.store $s/common.spher $s/common.print ; $f2 -c $opt $s/mofromao.f
prephss.o : $s/prephss.f $s/common.store ; $f2 -c $opt $s/prephss.f
prephdd.o : $s/prephdd.f $s/common.store ; $f2 -c $opt $s/prephdd.f
prtit.o : $s/prtit.f $s/common.store ; $f2 -c $opt $s/prtit.f
puttoiprop.o : $s/puttoiprop.f $s/common.basis $s/common.input ; $f2 -c $opt $s/puttoiprop.f
recoverwf.o : $s/recoverwf.f $s/common.store ; $f2 -c $opt $s/recoverwf.f
make1dm.o : $s/make1dm.f ; $f2 -c $opt $s/make1dm.f
make4xdm.o : $s/make4xdm.f ; $f2 -c $opt $s/make4xdm.f
makeguess.o : $s/makeguess.f $s/common.store $s/common.print ; $f2 -c $optno $s/makeguess.f
maptoold.o : $s/maptoold.f $s/common.store ; $f2 -c $opt $s/maptoold.f
mkdoub.o : $s/mkdoub.f $s/common.store ; $f2 -c $opt $s/mkdoub.f
mkdoubuc.o : $s/mkdoubuc.f $s/common.store ; $f2 -c $opt $s/mkdoubuc.f
mksing.o : $s/mksing.f $s/common.store ; $f2 -c $opt $s/mksing.f
ucdoubop.o : $s/ucdoubop.f $s/common.store ; $f2 -c $opt $s/ucdoubop.f
loaddoubb.o : $s/loaddoubb.f $s/common.store ; $f2 -c $opt $s/loaddoubb.f
loaddoubc.o : $s/loaddoubc.f ; $f2 -c $opt $s/loaddoubc.f
loadsinga.o : $s/loadsinga.f ; $f2 -c $opt $s/loadsinga.f
loadsingb.o : $s/loadsingb.f $s/common.store ; $f2 -c $opt $s/loadsingb.f
m12gen.o : $s/m12gen.f ; $f2 -c $opt $s/m12gen.f
m1gen.o : $s/m1gen.f ; $f2 -c $opt $s/m1gen.f
int1tobit.o : $s/int1tobit.f $s/common.store ; $f2 -c $opt $s/int1tobit.f
intcsf.o : $s/intcsf.f $s/common.basis $s/common.input $s/common.store $s/common.mrci $s/common.print ; $f2 -c $opt $s/intcsf.f
intcsfcas.o : $s/intcsfcas.f $s/common.basis $s/common.input $s/common.store $s/common.print ; $f2 -c $opt $s/intcsfcas.f
intsum.o : $s/intsum.f ; $f2 -c $opt $s/intsum.f
loaddouba.o : $s/loaddouba.f ; $f2 -c $opt $s/loaddouba.f
pslzzcsfcas.o : $s/pslzzcsfcas.f $s/common.store ; $f2 -c $opt $s/pslzzcsfcas.f
pythag.o : $s/pythag.f ; $f2 -c $opt $s/pythag.f
r2sb.o : $s/r2sb.f ; $f2 -c $opt $s/r2sb.f
strip2ps.o : $s/strip2ps.f $s/common.store ; $f2 -c $opt $s/strip2ps.f
stuffintohiicsf.o : $s/stuffintohiicsf.f $s/common.store ; $f2 -c $opt $s/stuffintohiicsf.f
updateiicsf.o : $s/updateiicsf.f $s/common.basis $s/common.input $s/common.store ; $f2 -c $opt $s/updateiicsf.f
updateq.o : $s/updateq.f ; $f2 -c $opt $s/updateq.f
lzziicsfcas.o : $s/lzziicsfcas.f $s/common.store ; $f2 -c $opt $s/lzziicsfcas.f
mov2ps.o : $s/mov2ps.f ; $f2 -c $opt $s/mov2ps.f
prodn.o : $s/prodn.f $s/common.store ; $f2 -c $opt $s/prodn.f
prtcsfvec.o : $s/prtcsfvec.f $s/common.mrci ; $f2 -c $opt $s/prtcsfvec.f
prtcsfveccas.o : $s/prtcsfveccas.f $s/common.store ; $f2 -c $opt $s/prtcsfveccas.f
pshamcsf.o : $s/pshamcsf.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/pshamcsf.f
pshamcsfcas.o : $s/pshamcsfcas.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/pshamcsfcas.f
pslzzcsf.o : $s/pslzzcsf.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/pslzzcsf.f
stuffintohii2.o : $s/stuffintohii2.f $s/common.store ; $f2 -c $opt $s/stuffintohii2.f
ilimts.o : $s/ilimts.f ; $f2 -c $opt $s/ilimts.f
int1tobit2.o : $s/int1tobit2.f $s/common.store ; $f2 -c $opt $s/int1tobit2.f
int1tobit3.o : $s/int1tobit3.f $s/common.store ; $f2 -c $opt $s/int1tobit3.f
lzziicsf.o : $s/lzziicsf.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/lzziicsf.f
uncompxu.o : $s/uncompxu.f ; $f2 -c $opt $s/uncompxu.f
hiicsf.o : $s/hiicsf.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hiicsf.f
hiicsfcas.o : $s/hiicsfcas.f $s/common.store $s/common.print ; $f2 -c $opt $s/hiicsfcas.f
hsingpart.o : $s/hsingpart.f $s/common.store ; $f2 -c $opt $s/hsingpart.f
hsscsf.o : $s/hsscsf.f $s/common.store $s/common.hf $s/common.mrci ; $f2 -c $opt $s/hsscsf.f
ifind2.o : $s/ifind2.f ; $f2 -c $opt $s/ifind2.f
hcis.o : $s/hcis.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hcis.f
hcsi.o : $s/hcsi.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hcsi.f
hcss.o : $s/hcss.f $s/common.store ; $f2 -c $opt $s/hcss.f
hddcsf.o : $s/hddcsf.f $s/common.store $s/common.hf $s/common.mrci ; $f2 -c $opt $s/hddcsf.f
hdoubpart.o : $s/hdoubpart.f $s/common.store ; $f2 -c $opt $s/hdoubpart.f
hdoubpartuc.o : $s/hdoubpartuc.f $s/common.store ; $f2 -c $opt $s/hdoubpartuc.f
hessycmo.o : $s/hessycmo.f $s/common.store ; $f2 -c $opt $s/hessycmo.f
hessyccas.o : $s/hessyccas.f $s/common.store ; $f2 -c $opt $s/hessyccas.f
hsscsf12.o : $s/hsscsf12.f $s/common.mrci $s/common.store ; $f2 -c $opt $s/hsscsf12.f
hcdduc.o : $s/hcdduc.f $s/common.store ; $f2 -c $opt $s/hcdduc.f
hcdi.o : $s/hcdi.f $s/common.store ; $f2 -c $opt $s/hcdi.f
hcdiuc.o : $s/hcdiuc.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hcdiuc.f
hcds.o : $s/hcds.f $s/common.store ; $f2 -c $opt $s/hcds.f
hcdsuc.o : $s/hcdsuc.f $s/common.store ; $f2 -c $optno $s/hcdsuc.f
hcid.o : $s/hcid.f $s/common.store ; $f2 -c $opt $s/hcid.f
hciduc.o : $s/hciduc.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hciduc.f
hddcsf12b.o : $s/hddcsf12b.f $s/common.mrci $s/common.store ; $f2 -c $opt $s/hddcsf12b.f
xtimesn.o : $s/xtimesn.f $s/common.store ; $f2 -c $opt $s/xtimesn.f
hccsfcas.o : $s/hccsfcas.f $s/common.store ; $f2 -c $opt $s/hccsfcas.f
hccsfd.o : $s/hccsfd.f $s/common.store $s/common.print ; $f2 -c $opt $s/hccsfd.f
hccsfrilr.o : $s/hccsfrilr.f $s/common.hf $s/common.store $s/common.mrci $s/common.print ; $f2 -c $opt $s/hccsfrilr.f
hcdd4v.o : $s/hcdd4v.f $s/common.store $s/common.mrci $s/common.hf ; $f2 -c $opt $s/hcdd4v.f
hcdddiag.o : $s/hcdddiag.f $s/common.store ; $f2 -c $opt $s/hcdddiag.f
hcddjk.o : $s/hcddjk.f $s/common.store ; $f2 -c $opt $s/hcddjk.f
prehcddjk.o : $s/prehcddjk.f $s/common.store ; $f2 -c $opt $s/prehcddjk.f
hcdduc1.o : $s/hcdduc1.f $s/common.store ; $f2 -c $opt $s/hcdduc1.f
hcdduc2.o : $s/hcdduc2.f $s/common.store ; $f2 -c $opt $s/hcdduc2.f
xtimesn2.o : $s/xtimesn2.f $s/common.store ; $f2 -c $opt $s/xtimesn2.f
hc1cd.o : $s/hc1cd.f $s/common.cas $s/common.store ; $f2 -c $opt $s/hc1cd.f
hc1cdbl.o : $s/hc1cdbl.f $s/common.cas $s/common.store ; $f2 -c $opt $s/hc1cdbl.f
hc1cnona1.o : $s/hc1cnona1.f $s/common.cas $s/common.store ; $f2 -c $opt $s/hc1cnona1.f
hc1cnona1d.o : $s/hc1cnona1d.f $s/common.cas $s/common.store ; $f2 -c $opt $s/hc1cnona1d.f
hccsf.o : $s/hccsf.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hccsf.f
mxmt.o : $s/mxmt.f ; $f2 -c $opt $s/mxmt.f
paraerik.o : $s/paraerik.f $s/common.hf $s/common.store $s/common.basis $s/common.mrci ; $f2 -c $opt $s/paraerik.f
sandtd.o : $s/sandtd.f $s/common.store ; $f2 -c $opt $s/sandtd.f
vdmo2so.o : $s/vdmo2so.f $s/common.store $s/common.hf $s/common.mrci ; $f2 -c $opt $s/vdmo2so.f
getbas0.o : $s/getbas0.f ; $f2 -c $opt $s/getbas0.f
getwf.o : $s/getwf.f $s/common.store $s/common.basis ; $f2 -c $optno $s/getwf.f
gotoqs.o : $s/gotoqs.f ; $f2 -c $opt $s/gotoqs.f
graborb.o : $s/graborb.f $s/common.store $s/common.basis ; $f2 -c $opt $s/graborb.f
grado.o : $s/grado.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $optno $s/grado.f
gradcas.o : $s/gradcas.f $s/common.hf $s/common.store $s/common.spher $s/common.basis $s/common.print ; $f2 -c $opt $s/gradcas.f
gtmom.o : $s/gtmom.f $s/common.store $s/common.basis ; $f2 -c $opt $s/gtmom.f
gtmomso.o : $s/gtmomso.f $s/common.store ; $f2 -c $opt $s/gtmomso.f
hc1c.o : $s/hc1c.f $s/common.cas $s/common.store ; $f2 -c $opt $s/hc1c.f
parajkfromhd2c.o : $s/parajkfromhd2c.f $s/common.store $s/common.hf ; $f2 -c $opt $s/parajkfromhd2c.f
psioppsi.o : $s/psioppsi.f $s/common.store ; $f2 -c $opt $s/psioppsi.f
rordr.o : $s/rordr.f $s/common.hf $s/common.store ; $f2 -c $opt $s/rordr.f
sym4o.o : $s/sym4o.f $s/common.store ; $f2 -c $opt $s/sym4o.f
tofrob.o : $s/tofrob.f $s/common.store ; $f2 -c $opt $s/tofrob.f
transder.o : $s/transder.f $s/common.store $s/common.hf ; $f2 -c $opt $s/transder.f
transder1.o : $s/transder1.f $s/common.store $s/common.hf ; $f2 -c $opt $s/transder1.f
trianglez.o : $s/trianglez.f $s/common.store ; $f2 -c $opt $s/trianglez.f
parajkfromhd1j.o : $s/parajkfromhd1j.f $s/common.store $s/common.hf ; $f2 -c $opt $s/parajkfromhd1j.f
parajkfromhd2b.o : $s/parajkfromhd2b.f $s/common.store $s/common.hf ; $f2 -c $opt $s/parajkfromhd2b.f
testme.o : $s/testme.f $s/common.store ; $f2 -c $opt $s/testme.f
unmakegd.o : $s/unmakegd.f $s/common.store ; $f2 -c $opt $s/unmakegd.f
unmakegs.o : $s/unmakegs.f $s/common.store ; $f2 -c $opt $s/unmakegs.f
parah0d2.o : $s/parah0d2.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parah0d2.f
parah0grad.o : $s/parah0grad.f $s/common.hf $s/common.store $s/common.spher $s/common.basis ; $f2 -c $opt $s/parah0grad.f
parajkfromhd0.o : $s/parajkfromhd0.f $s/common.store $s/common.hf ; $f2 -c $opt $s/parajkfromhd0.f
parajkfromhd1.o : $s/parajkfromhd1.f $s/common.store $s/common.hf ; $f2 -c $opt $s/parajkfromhd1.f
sotest2.o : $s/sotest2.f $s/common.store ; $f2 -c $opt $s/sotest2.f
sotest2c.o : $s/sotest2c.f $s/common.store ; $f2 -c $opt $s/sotest2c.f
sotest.o : $s/sotest.f $s/common.store ; $f2 -c $opt $s/sotest.f
onedints.o : $s/onedints.f $s/common.store ; $f2 -c $opt $s/onedints.f
orbder.o : $s/orbder.f $s/common.store $s/common.hf ; $f2 -c $opt $s/orbder.f
orbdercas.o : $s/orbdercas.f $s/common.store $s/common.hf ; $f2 -c $opt $s/orbdercas.f
paraerid2b.o : $s/paraerid2b.f $s/common.hf $s/common.store $s/common.basis ; $f2 -c $opt $s/paraerid2b.f
paraerid2c.o : $s/paraerid2c.f $s/common.hf $s/common.store $s/common.basis ; $f2 -c $opt $s/paraerid2c.f
paraerid.o : $s/paraerid.f $s/common.hf $s/common.store $s/common.basis ; $f2 -c $opt $s/paraerid.f
paraeridj.o : $s/paraeridj.f $s/common.hf $s/common.store $s/common.basis ; $f2 -c $opt $s/paraeridj.f
int4copy.o : $s/int4copy.f $s/common.store ; $f2 -c $opt $s/int4copy.f
makegd.o : $s/makegd.f ; $f2 -c $opt $s/makegd.f
makegs.o : $s/makegs.f $s/common.store ; $f2 -c $opt $s/makegs.f
mtimesx.o : $s/mtimesx.f $s/common.store ; $f2 -c $opt $s/mtimesx.f
nattracd.o : $s/nattracd.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/nattracd.f
oned2ints.o : $s/oned2ints.f $s/common.store ; $f2 -c $opt $s/oned2ints.f
oneider.o : $s/oneider.f $s/common.store $s/common.spher ; $f2 -c $opt $s/oneider.f
restid.o : $s/restid.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/restid.f
restind.o : $s/restind.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/restind.f
genwf.o : $s/genwf.f $s/common.store ; $f2 -c $optno $s/genwf.f
hcbk.o : $s/hcbk.f $s/common.store ; $f2 -c $opt $s/hcbk.f
hccsfbk.o : $s/hccsfbk.f $s/common.store $s/common.print ; $f2 -c $opt $s/hccsfbk.f
hccsfd1.o : $s/hccsfd1.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hccsfd1.f
hcddjkd1.o : $s/hcddjkd1.f $s/common.store ; $f2 -c $opt $s/hcddjkd1.f
hcisd1.o : $s/hcisd1.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hcisd1.f
hcssd1.o : $s/hcssd1.f $s/common.store ; $f2 -c $opt $s/hcssd1.f
hcdsbk.o : $s/hcdsbk.f $s/common.store ; $f2 -c $opt $s/hcdsbk.f
hcidbk.o : $s/hcidbk.f $s/common.store ; $f2 -c $opt $s/hcidbk.f
hcidbkuc.o : $s/hcidbkuc.f $s/common.store ; $f2 -c $opt $s/hcidbkuc.f
hcisbk.o : $s/hcisbk.f $s/common.store ; $f2 -c $opt $s/hcisbk.f
hcsdbk.o : $s/hcsdbk.f $s/common.store ; $f2 -c $opt $s/hcsdbk.f
hcsibk.o : $s/hcsibk.f $s/common.store ; $f2 -c $opt $s/hcsibk.f
hcssbk.o : $s/hcssbk.f $s/common.store ; $f2 -c $opt $s/hcssbk.f
wnorm.o : $s/wnorm.f $s/common.store ; $f2 -c $opt $s/wnorm.f
get4int.o : $s/get4int.f $s/common.store ; $f2 -c $opt $s/get4int.f
hcddjkbk.o : $s/hcddjkbk.f $s/common.store ; $f2 -c $opt $s/hcddjkbk.f
hcdducbk.o : $s/hcdducbk.f $s/common.store ; $f2 -c $opt $s/hcdducbk.f
hcdibk.o : $s/hcdibk.f $s/common.store ; $f2 -c $opt $s/hcdibk.f
hcdibkuc.o : $s/hcdibkuc.f $s/common.store ; $f2 -c $opt $s/hcdibkuc.f
hcdsd1.o : $s/hcdsd1.f $s/common.store ; $f2 -c $opt $s/hcdsd1.f
wnorm1.o : $s/wnorm1.f $s/common.store ; $f2 -c $opt $s/wnorm1.f
wnorm2.o : $s/wnorm2.f $s/common.store ; $f2 -c $opt $s/wnorm2.f
wnorm2c.o : $s/wnorm2c.f ; $f2 -c $opt $s/wnorm2c.f
genryd.o : $s/genryd.f $s/common.hf $s/common.store $s/common.cas $s/common.print ; $f2 -c $opt $s/genryd.f
genmatn.o : $s/genmatn.f $s/common.store ; $f2 -c $opt $s/genmatn.f
genmatn2.o : $s/genmatn2.f $s/common.store ; $f2 -c $opt $s/genmatn2.f
genmatn3.o : $s/genmatn3.f $s/common.store ; $f2 -c $opt $s/genmatn3.f
singnm1.o : $s/singnm1.f $s/common.store ; $f2 -c $opt $s/singnm1.f
genmat.o : $s/genmat.f $s/common.store ; $f2 -c $opt $s/genmat.f
genhsoa2.o : $s/genhsoa2.f $s/common.store ; $f2 -c $opt $s/genhsoa2.f
genfcnp.o : $s/genfcnp.f $s/common.store $s/common.mrci $s/common.print ; $f2 -c $opt $s/genfcnp.f
psisopsi.o : $s/psisopsi.f $s/common.store $s/common.print ; $f2 -c $opt $s/psisopsi.f
genfcn2.o : $s/genfcn2.f $s/common.store $s/common.mrci $s/common.print ; $f2 -c $opt $s/genfcn2.f
genf0.o : $s/genf0.f $s/common.store ; $f2 -c $opt $s/genf0.f
genergy.o : $s/genergy.f $s/common.store $s/common.print ; $f2 -c $opt $s/genergy.f
gendertype.o : $s/gendertype.f $s/common.hf ; $f2 -c $opt $s/gendertype.f
gencsf2b.o : $s/gencsf2b.f $s/common.store $s/common.print ; $f2 -c $opt $s/gencsf2b.f
gencsf3.o : $s/gencsf3.f $s/common.store $s/common.print ; $f2 -c $opt $s/gencsf3.f
hcidbk4.o : $s/hcidbk4.f $s/common.store ; $f2 -c $opt $s/hcidbk4.f
hcisbk4.o : $s/hcisbk4.f $s/common.store ; $f2 -c $opt $s/hcisbk4.f
hcsdbk4.o : $s/hcsdbk4.f $s/common.store ; $f2 -c $opt $s/hcsdbk4.f
hcsibk4.o : $s/hcsibk4.f $s/common.store ; $f2 -c $opt $s/hcsibk4.f
hcssbk4.o : $s/hcssbk4.f $s/common.store ; $f2 -c $opt $s/hcssbk4.f
i4to8copy.o : $s/i4to8copy.f ; $f2 -c $opt $s/i4to8copy.f
paraerikbk4.o : $s/paraerikbk4.f $s/common.hf $s/common.store $s/common.basis $s/common.mrci ; $f2 -c $opt $s/paraerikbk4.f
testsome.o : $s/testsome.f $s/common.store ; $f2 -c $opt $s/testsome.f
vdmo2sobk4.o : $s/vdmo2sobk4.f $s/common.store $s/common.hf $s/common.mrci ; $f2 -c $opt $s/vdmo2sobk4.f
getcup.o : $s/getcup.f $s/common.store ; $f2 -c $opt $s/getcup.f
hccsfbk4.o : $s/hccsfbk4.f $s/common.store $s/common.print ; $f2 -c $opt $s/hccsfbk4.f
hcddjkbk4.o : $s/hcddjkbk4.f $s/common.store ; $f2 -c $opt $s/hcddjkbk4.f
hcdibk4.o : $s/hcdibk4.f $s/common.store ; $f2 -c $opt $s/hcdibk4.f
hcdsbk4.o : $s/hcdsbk4.f $s/common.store ; $f2 -c $opt $s/hcdsbk4.f
spinloop.o : $s/spinloop.f $s/common.store ; $f2 -c $opt $s/spinloop.f
spinloop1.o : $s/spinloop1.f $s/common.store ; $f2 -c $opt $s/spinloop1.f
gen12nona1.o : $s/gen12nona1.f $s/common.store ; $f2 -c $opt $s/gen12nona1.f
gencsf1.o : $s/gencsf1.f ; $f2 -c $opt $s/gencsf1.f
gencsf2.o : $s/gencsf2.f $s/common.store $s/common.print ; $f2 -c $opt $s/gencsf2.f
gencup.o : $s/gencup.f $s/common.store ; $f2 -c $opt $s/gencup.f
uncompxut.o : $s/uncompxut.f ; $f2 -c $opt $s/uncompxut.f
gcode.o : $s/gcode.f ; $f2 -c $opt $s/gcode.f
gcsfps.o : $s/gcsfps.f ; $f2 -c $opt $s/gcsfps.f
gcsfpslz.o : $s/gcsfpslz.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/gcsfpslz.f
gen12.o : $s/gen12.f $s/common.store ; $f2 -c $opt $s/gen12.f
epdvec.o : $s/epdvec.f $s/common.store ; $f2 -c $opt $s/epdvec.f
erir.o : $s/erir.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/erir.f
erird.o : $s/erird.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/erird.f
fiddleh.o : $s/fiddleh.f $s/common.store ; $f2 -c $opt $s/fiddleh.f
fillrup.o : $s/fillrup.f $s/common.store ; $f2 -c $opt $s/fillrup.f
gandc.o : $s/gandc.f $s/common.store ; $f2 -c $opt $s/gandc.f
gandcr.o : $s/gandcr.f $s/common.store ; $f2 -c $opt $s/gandcr.f
gandc4.o : $s/gandc4.f $s/common.store ; $f2 -c $opt $s/gandc4.f
dynwtr.o : $s/dynwtr.f ; $f2 -c $opt $s/dynwtr.f
dumpbas0.o : $s/dumpbas0.f ; $f2 -c $opt $s/dumpbas0.f
dumpbas2.o : $s/dumpbas2.f ; $f2 -c $opt $s/dumpbas2.f
dumpbas.o : $s/dumpbas.f $s/common.store ; $f2 -c $opt $s/dumpbas.f
dumpdoubc.o : $s/dumpdoubc.f $s/common.store ; $f2 -c $opt $s/dumpdoubc.f
dumpdoubuc.o : $s/dumpdoubuc.f $s/common.store ; $f2 -c $opt $s/dumpdoubuc.f
dumph.o : $s/dumph.f ; $f2 -c $opt $s/dumph.f
dumpsing.o : $s/dumpsing.f $s/common.store ; $f2 -c $opt $s/dumpsing.f
dosingnona1.o : $s/dosingnona1.f $s/common.store ; $f2 -c $opt $s/dosingnona1.f
dosingnona1d.o : $s/dosingnona1d.f $s/common.store ; $f2 -c $opt $s/dosingnona1d.f
dotdvecs.o : $s/dotdvecs.f $s/common.print ; $f2 -c $opt $s/dotdvecs.f
dotvall.o : $s/dotvall.f $s/common.store ; $f2 -c $opt $s/dotvall.f
dovr.o : $s/dovr.f $s/common.store ; $f2 -c $opt $s/dovr.f
dumb.o : $s/dumb.f ; $f2 -c $opt $s/dumb.f
dumpbas1.o : $s/dumpbas1.f ; $f2 -c $opt $s/dumpbas1.f
doab2dbl.o : $s/doab2dbl.f $s/common.store ; $f2 -c $opt $s/doab2dbl.f
dodoub.o : $s/dodoub.f $s/common.store ; $f2 -c $opt $s/dodoub.f
dodoubd.o : $s/dodoubd.f $s/common.store ; $f2 -c $opt $s/dodoubd.f
dodoubdbl.o : $s/dodoubdbl.f $s/common.store ; $f2 -c $opt $s/dodoubdbl.f
dodoubnona1.o : $s/dodoubnona1.f $s/common.store ; $f2 -c $opt $s/dodoubnona1.f
dodoubnona1d.o : $s/dodoubnona1d.f $s/common.store ; $f2 -c $opt $s/dodoubnona1d.f
dosing.o : $s/dosing.f $s/common.store ; $f2 -c $opt $s/dosing.f
dosingd.o : $s/dosingd.f $s/common.store ; $f2 -c $opt $s/dosingd.f
dosingdbl.o : $s/dosingdbl.f $s/common.store ; $f2 -c $opt $s/dosingdbl.f
doab11.o : $s/doab11.f $s/common.store ; $f2 -c $opt $s/doab11.f
doab1.o : $s/doab1.f $s/common.store ; $f2 -c $opt $s/doab1.f
doab1d.o : $s/doab1d.f $s/common.store ; $f2 -c $opt $s/doab1d.f
doab1dbl.o : $s/doab1dbl.f $s/common.store ; $f2 -c $opt $s/doab1dbl.f
doab1non.o : $s/doab1non.f $s/common.store ; $f2 -c $opt $s/doab1non.f
doab1nond.o : $s/doab1nond.f $s/common.store ; $f2 -c $opt $s/doab1nond.f
doab2.o : $s/doab2.f $s/common.store ; $f2 -c $opt $s/doab2.f
doab2d.o : $s/doab2d.f $s/common.store ; $f2 -c $opt $s/doab2d.f
doab21.o : $s/doab21.f $s/common.store ; $f2 -c $opt $s/doab21.f
diagx.o : $s/diagx.f $s/common.store $s/common.hf ; $f2 -c $opt $s/diagx.f
diagy.o : $s/diagy.f $s/common.store $s/common.hf ; $f2 -c $opt $s/diagy.f
distit8.o : $s/distit8.f ; $f2 -c $opt $s/distit8.f
distit.o : $s/distit.f ; $f2 -c $opt $s/distit.f
do4v.o : $s/do4v.f $s/common.store ; $f2 -c $opt $s/do4v.f
derofxor.o : $s/derofxor.f $s/common.hf $s/common.store $s/common.spher ; $f2 -c $opt $s/derofxor.f
dgerid.o : $s/dgerid.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/dgerid.f
diagdy.o : $s/diagdy.f $s/common.store $s/common.print ; $f2 -c $opt $s/diagdy.f
diaghii.o : $s/diaghii.f $s/common.store $s/common.cas ; $f2 -c $opt $s/diaghii.f
diaghiid.o : $s/diaghiid.f $s/common.store $s/common.cas ; $f2 -c $opt $s/diaghiid.f
diaghiidbl.o : $s/diaghiidbl.f $s/common.store $s/common.cas ; $f2 -c $opt $s/diaghiidbl.f
diaglzz.o : $s/diaglzz.f $s/common.store $s/common.cas ; $f2 -c $opt $s/diaglzz.f
derh01b.o : $s/derh01b.f $s/common.store ; $f2 -c $opt $s/derh01b.f
derh0.o : $s/derh0.f $s/common.hf $s/common.store $s/common.spher ; $f2 -c $opt $s/derh0.f
derid.o : $s/derid.f $s/common.store $s/common.spher $s/common.rys ; $f2 -c $opt $s/derid.f
derofxor1.o : $s/derofxor1.f $s/common.hf $s/common.store $s/common.spher ; $f2 -c $opt $s/derofxor1.f
diaghii1.o : $s/diaghii1.f $s/common.store $s/common.cas $s/common.hf ; $f2 -c $opt $s/diaghii1.f
diaglzz1.o : $s/diaglzz1.f $s/common.store $s/common.cas ; $f2 -c $opt $s/diaglzz1.f
deconvert2torc.o : $s/deconvert2torc.f $s/common.store ; $f2 -c $opt $s/deconvert2torc.f
deconvert12torc.o : $s/deconvert12torc.f $s/common.store ; $f2 -c $opt $s/deconvert12torc.f
denmx.o : $s/denmx.f $s/common.store ; $f2 -c $opt $s/denmx.f
denmx12.o : $s/denmx12.f $s/common.hf $s/common.store ; $f2 -c $opt $s/denmx12.f
der3part.o : $s/der3part.f $s/common.store $s/common.hf ; $f2 -c $opt $s/der3part.f
deramat.o : $s/deramat.f $s/common.hf $s/common.store ; $f2 -c $opt $s/deramat.f
derh01.o : $s/derh01.f $s/common.basis $s/common.hf $s/common.store $s/common.spher ; $f2 -c $opt $s/derh01.f
hccsfd12.o : $s/hccsfd12.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hccsfd12.f
hcsid12.o : $s/hcsid12.f $s/common.store $s/common.mrci ; $f2 -c $opt $s/hcsid12.f
corenota1.o : $s/corenota1.f $s/common.store $s/common.hf $s/common.cas $s/common.print ; $f2 -c $opt $s/corenota1.f
corelz2.o : $s/corelz2.f $s/common.store $s/common.hf $s/common.cas ; $f2 -c $opt $s/corelz2.f
countcc.o : $s/countcc.f ; $f2 -c $opt $s/countcc.f
countem2.o : $s/countem2.f $s/common.mrci $s/common.store ; $f2 -c $opt $s/countem2.f
countfcn.o : $s/countfcn.f ; $f2 -c $opt $s/countfcn.f
cupo1rp.o : $s/cupo1rp.f $s/common.store $s/common.print ; $f2 -c $opt $s/cupo1rp.f
cupo21.o : $s/cupo21.f $s/common.store $s/common.print ; $f2 -c $opt $s/cupo21.f
dcannon.o : $s/dcannon.f $s/common.hf $s/common.store $s/common.print ; $f2 -c $opt $s/dcannon.f
dcbit.o : $s/dcbit.f ; $f2 -c $opt $s/dcbit.f
cont2ecsf.o : $s/cont2ecsf.f $s/common.mrci $s/common.store $s/common.print ; $f2 -c $opt $s/cont2ecsf.f
cont2ecsfb.o : $s/cont2ecsfb.f $s/common.store $s/common.basis $s/common.input $s/common.mrci $s/common.print ; $f2 -c $opt $s/cont2ecsfb.f
cont2ecsfa.o : $s/cont2ecsfa.f $s/common.store $s/common.mrci $s/common.print ; $f2 -c $opt $s/cont2ecsfa.f
convert12tor.o : $s/convert12tor.f $s/common.store ; $f2 -c $opt $s/convert12tor.f
convert12torc.o : $s/convert12torc.f $s/common.store ; $f2 -c $opt $s/convert12torc.f
convert1tor.o : $s/convert1tor.f $s/common.store ; $f2 -c $opt $s/convert1tor.f
convert2tor.o : $s/convert2tor.f $s/common.store ; $f2 -c $opt $s/convert2tor.f
convert2torc.o : $s/convert2torc.f $s/common.store ; $f2 -c $opt $s/convert2torc.f
core.o : $s/core.f $s/common.store $s/common.hf $s/common.cas $s/common.print ; $f2 -c $opt $s/core.f
consolidate1.o : $s/consolidate1.f $s/common.store ; $f2 -c $opt $s/consolidate1.f
matchc.o : $s/matchc.f $s/common.store $s/common.print ; $f2 -c $opt $s/matchc.f
orthovrb.o : $s/orthovrb.f $s/common.store $s/common.print ; $f2 -c $opt $s/orthovrb.f
compxtimes.o : $s/compxtimes.f ; $f2 -c $opt $s/compxtimes.f
consolidate3.o : $s/consolidate3.f $s/common.store ; $f2 -c $opt $s/consolidate3.f
consolidate2.o : $s/consolidate2.f ; $f2 -c $opt $s/consolidate2.f
matchcl.o : $s/matchcl.f $s/common.print ; $f2 -c $opt $s/matchcl.f
check2.o : $s/check2.f $s/common.store ; $f2 -c $opt $s/check2.f
check1.o : $s/check1.f ; $f2 -c $opt $s/check1.f
check3.o : $s/check3.f ; $f2 -c $opt $s/check3.f
checkc.o : $s/checkc.f $s/common.store ; $f2 -c $opt $s/checkc.f
checkps.o : $s/checkps.f $s/common.store ; $f2 -c $opt $s/checkps.f
cmpvcsf.o : $s/cmpvcsf.f ; $f2 -c $opt $s/cmpvcsf.f
cmpvcsf2.o : $s/cmpvcsf2.f ; $f2 -c $opt $s/cmpvcsf2.f
cnt12.o : $s/cnt12.f ; $f2 -c $opt $s/cnt12.f
cnt12nona1.o : $s/cnt12nona1.f ; $f2 -c $opt $s/cnt12nona1.f
compxtimes2.o : $s/compxtimes2.f $s/common.store ; $f2 -c $opt $s/compxtimes2.f
bodcpart.o : $s/bodcpart.f $s/common.store ; $f2 -c $opt $s/bodcpart.f
buildder3s.o : $s/buildder3s.f $s/common.store $s/common.hf ; $f2 -c $opt $s/buildder3s.f
buildhess.o : $s/buildhess.f $s/common.hf $s/common.store ; $f2 -c $opt $s/buildhess.f
buildhessdx.o : $s/buildhessdx.f $s/common.hf $s/common.store ; $f2 -c $opt $s/buildhessdx.f
c2eprop.o : $s/c2eprop.f $s/common.store $s/common.basis ; $f2 -c $opt $s/c2eprop.f
cal1int.o : $s/cal1int.f $s/common.hf $s/common.store ; $f2 -c $opt $s/cal1int.f
cart2spher.o : $s/cart2spher.f $s/common.spher $s/common.store ; $f2 -c $opt $s/cart2spher.f
cfileget.o : $s/cfileget.f $s/common.store ; $f2 -c $opt $s/cfileget.f
atomll.o : $s/atomll.f $s/common.store ; $f2 -c $opt $s/atomll.f
avgr.o : $s/avgr.f ; $f2 -c $opt $s/avgr.f
mapv2v.o : $s/mapv2v.f $s/common.store ; $f2 -c $opt $s/mapv2v.f
pdump.o : $s/pdump.f ; $f2 -c $opt $s/pdump.f
 addcomma8.o : $s/addcomma8.f ; $f2 -c $opt $s/addcomma8.f
 addqtop.o : $s/addqtop.f $s/common.store ; $f2 -c $opt $s/addqtop.f
!
echo "   59 format(/'Linked on: " `date` "')" > $s/mpec.date
echo "      pwd=" >> $s/mpec.date
echo "     $'" `pwd` "'" >> $s/mpec.date
sleep 5
make -f makefile mpec
chmod og+rx mpec
mv mpec ..

