# Results makefile

ILIB = $(LIBDIR)/$(FLIB)
BINDIR = ../$(FBIN)

all :	xairy.reslt xamebsa.reslt xamoeba.reslt xanneal.reslt \
	xarcode.reslt xavevar.reslt xbalanc.reslt xbandec.reslt \
	xbanmul.reslt xbcucof.reslt xbcuint.reslt xbeschb.reslt \
	xbessi.reslt xbessi0.reslt xbessi1.reslt xbessik.reslt \
	xbessj.reslt xbessj0.reslt xbessj1.reslt xbessjy.reslt \
	xbessk.reslt xbessk0.reslt xbessk1.reslt xbessy.reslt \
	xbessy0.reslt xbessy1.reslt xbeta.reslt xbetai.reslt xbico.reslt \
	xbnldev.reslt xbrent.reslt xbroydn.reslt xbsstep.reslt \
	xcaldat.reslt xchder.reslt xchebev.reslt xchebft.reslt \
	xchebpc.reslt xchint.reslt xcholsl.reslt xchsone.reslt \
	xchstwo.reslt xcisi.reslt xcntab1.reslt xcntab2.reslt \
	xconvlv.reslt xcorrel.reslt xcosft1.reslt xcosft2.reslt \
	xcovsrt.reslt xcrank.reslt xcyclic.reslt xdawson.reslt \
	xdbrent.reslt xddpoly.reslt xdecchk.reslt xdf1dim.reslt \
	xdfpmin.reslt xdfridr.reslt xdftint.reslt xeclass.reslt \
	xeclazz.reslt xei.reslt xeigsrt.reslt xelle.reslt xellf.reslt \
	xellpi.reslt xelmhes.reslt xerf.reslt xerfc.reslt xerfcc.reslt \
	xeulsum.reslt xevlmem.reslt xexpdev.reslt xexpint.reslt \
	xf1dim.reslt xfactln.reslt xfactrl.reslt xfasper.reslt \
	xfgauss.reslt xfit.reslt xfitexy.reslt xfixrts.reslt xfleg.reslt \
	xflmoon.reslt xfour1.reslt xfourfs.reslt xfourn.reslt \
	xfpoly.reslt xfred2.reslt xfredin.reslt xfrenel.reslt \
	xfrprmn.reslt xftest.reslt xgamdev.reslt xgammln.reslt \
	xgammp.reslt xgammq.reslt xgasdev.reslt xgaucof.reslt \
	xgauher.reslt xgaujac.reslt xgaulag.reslt xgauleg.reslt \
	xgaussj.reslt xgcf.reslt xgolden.reslt xgser.reslt \
	xhpsel.reslt xhpsort.reslt xhqr.reslt xhuffman.reslt xhunt.reslt \
	xhypgeo.reslt xicrc.reslt xigray.reslt xindexx.reslt \
	xirbit1.reslt xirbit2.reslt xjacobi.reslt xjulday.reslt \
	xkendl1.reslt xkendl2.reslt xks2d1s.reslt xks2d2s.reslt \
	xksone.reslt xkstwo.reslt xlaguer.reslt xlfit.reslt \
	xlinbcg.reslt xlinmin.reslt xlocate.reslt xlubksb.reslt \
	xludcmp.reslt xmachar.reslt xmedfit.reslt xmemcof.reslt \
	xmgfas.reslt xmglin.reslt xmidpnt.reslt xmiser.reslt xmmid.reslt \
	xmnbrak.reslt xmnewt.reslt xmoment.reslt xmppi.reslt \
	xmprove.reslt xmrqcof.reslt xmrqmin.reslt xnewt.reslt \
	xodeint.reslt xorthog.reslt xpade.reslt xpccheb.reslt \
	xpcshft.reslt xpearsn.reslt xperiod.reslt xpiksr2.reslt \
	xpiksrt.reslt xplgndr.reslt xpoidev.reslt xpolcoe.reslt \
	xpolcof.reslt xpoldiv.reslt xpolin2.reslt xpolint.reslt \
	xpowell.reslt xpredic.reslt xprobks.reslt xpsdes.reslt \
	xpzextr.reslt xqgaus.reslt xqrdcmp.reslt xqromb.reslt \
	xqromo.reslt xqroot.reslt xqrsolv.reslt xqrupdt.reslt \
	xqsimp.reslt xqtrap.reslt xquad3d.reslt xran.reslt xran4.reslt \
	xrank.reslt xratint.reslt xratlsq.reslt xrc.reslt xrd.reslt \
	xrealft.reslt xrf.reslt xrj.reslt xrk4.reslt xrkdumb.reslt \
	xrkqs.reslt xrlft3.reslt xrofunc.reslt xrtbis.reslt \
	xrtflsp.reslt xrtnewt.reslt xrtsafe.reslt xrtsec.reslt \
	xrzextr.reslt xsavgol.reslt xscrsho.reslt xselect.reslt \
	xselip.reslt xshell.reslt xsimplx.reslt xsimpr.reslt \
	xsinft.reslt xsncndn.reslt xsobseq.reslt xsor.reslt xsort.reslt \
	xsort2.reslt xsort3.reslt xspctrm.reslt xspear.reslt \
	xsphbes.reslt xsphfpt.reslt xsplie2.reslt xsplin2.reslt \
	xspline.reslt xsplint.reslt xsprsax.reslt xsprsin.reslt \
	xsprspm.reslt xsprstm.reslt xsprstp.reslt xsprstx.reslt \
	xstifbs.reslt xstiff.reslt xstoerm.reslt xsvbksb.reslt \
	xsvdcmp.reslt xsvdfit.reslt xsvdvar.reslt xtoeplz.reslt \
	xtptest.reslt xtqli.reslt xtrapzd.reslt xtred2.reslt \
	xtridag.reslt xttest.reslt xtutest.reslt xtwofft.reslt \
	xvander.reslt xvegas.reslt xvoltra.reslt xwt1.reslt xwtn.reslt \
	xzbrac.reslt xzbrak.reslt xzbrent.reslt xzrhqr.reslt \
	xzriddr.reslt xzroots.reslt badluk.reslt fredex.reslt \
	sfroid.reslt sphoot.reslt

# The following is extremely inelegant and ugly, but there seems to be no
# better way to accomplish the objectives in a manner portable across
# a wide variety of systems

xairy.reslt : ../src/xairy.f $(ILIB)
	-if [ -r ../bin/xairy.out ]; then \
	../bin/xairy.out > xairy.reslt; \
	else \
	$(FC) $(FFLAGS) -o xairy.out ../src/xairy.f $(ILIB) && \
	xairy.out > xairy.reslt && /bin/rm -f xairy.o xairy.out; \
	fi

xamebsa.reslt : ../src/xamebsa.f $(ILIB)
	-if [ -r ../bin/xamebsa.out ]; then \
	../bin/xamebsa.out < ../responses/xamebsa.resp > xamebsa.reslt; \
	else \
	$(FC) $(FFLAGS) -o xamebsa.out ../src/xamebsa.f $(ILIB) && \
	xamebsa.out < ../responses/xamebsa.resp > xamebsa.reslt && \
	/bin/rm -f xamebsa.o xamebsa.out; \
	fi

xamoeba.reslt : ../src/xamoeba.f $(ILIB)
	-if [ -r ../bin/xamoeba.out ]; then \
	../bin/xamoeba.out > xamoeba.reslt; \
	else \
	$(FC) $(FFLAGS) -o xamoeba.out ../src/xamoeba.f $(ILIB) && \
	xamoeba.out > xamoeba.reslt && /bin/rm -f xamoeba.o xamoeba.out; \
	fi

xanneal.reslt : ../src/xanneal.f $(ILIB)
	-if [ -r ../bin/xanneal.out ]; then \
	../bin/xanneal.out > xanneal.reslt; \
	else \
	$(FC) $(FFLAGS) -o xanneal.out ../src/xanneal.f $(ILIB) && \
	xanneal.out > xanneal.reslt && /bin/rm -f xanneal.o xanneal.out; \
	fi

xarcode.reslt : ../src/xarcode.f $(ILIB)
	-if [ -r ../bin/xarcode.out ]; then \
	../bin/xarcode.out < ../responses/xarcode.resp > xarcode.reslt; \
	else \
	$(FC) $(FFLAGS) -o xarcode.out ../src/xarcode.f $(ILIB) && \
	xarcode.out < ../responses/xarcode.resp > xarcode.reslt && \
	/bin/rm -f xarcode.o xarcode.out; \
	fi

xavevar.reslt : ../src/xavevar.f $(ILIB)
	-if [ -r ../bin/xavevar.out ]; then \
	../bin/xavevar.out > xavevar.reslt; \
	else \
	$(FC) $(FFLAGS) -o xavevar.out ../src/xavevar.f $(ILIB) && \
	xavevar.out > xavevar.reslt && /bin/rm -f xavevar.o xavevar.out; \
	fi

xbalanc.reslt : ../src/xbalanc.f $(ILIB)
	-if [ -r ../bin/xbalanc.out ]; then \
	../bin/xbalanc.out > xbalanc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbalanc.out ../src/xbalanc.f $(ILIB) && \
	xbalanc.out > xbalanc.reslt && /bin/rm -f xbalanc.o xbalanc.out; \
	fi

xbandec.reslt : ../src/xbandec.f $(ILIB)
	-if [ -r ../bin/xbandec.out ]; then \
	../bin/xbandec.out > xbandec.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbandec.out ../src/xbandec.f $(ILIB) && \
	xbandec.out > xbandec.reslt && /bin/rm -f xbandec.o xbandec.out; \
	fi

xbanmul.reslt : ../src/xbanmul.f $(ILIB)
	-if [ -r ../bin/xbanmul.out ]; then \
	../bin/xbanmul.out > xbanmul.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbanmul.out ../src/xbanmul.f $(ILIB) && \
	xbanmul.out > xbanmul.reslt && /bin/rm -f xbanmul.o xbanmul.out; \
	fi

xbcucof.reslt : ../src/xbcucof.f $(ILIB)
	-if [ -r ../bin/xbcucof.out ]; then \
	../bin/xbcucof.out > xbcucof.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbcucof.out ../src/xbcucof.f $(ILIB) && \
	xbcucof.out > xbcucof.reslt && /bin/rm -f xbcucof.o xbcucof.out; \
	fi

xbcuint.reslt : ../src/xbcuint.f $(ILIB)
	-if [ -r ../bin/xbcuint.out ]; then \
	../bin/xbcuint.out > xbcuint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbcuint.out ../src/xbcuint.f $(ILIB) && \
	xbcuint.out > xbcuint.reslt && /bin/rm -f xbcuint.o xbcuint.out; \
	fi

xbeschb.reslt : ../src/xbeschb.f $(ILIB)
	-if [ -r ../bin/xbeschb.out ]; then \
	../bin/xbeschb.out < ../responses/xbeschb.resp > xbeschb.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbeschb.out ../src/xbeschb.f $(ILIB) && \
	xbeschb.out < ../responses/xbeschb.resp > xbeschb.reslt && \
	/bin/rm -f xbeschb.o xbeschb.out; \
	fi

xbessi.reslt : ../src/xbessi.f $(ILIB)
	-if [ -r ../bin/xbessi.out ]; then \
	../bin/xbessi.out > xbessi.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessi.out ../src/xbessi.f $(ILIB) && \
	xbessi.out > xbessi.reslt && /bin/rm -f xbessi.o xbessi.out; \
	fi

xbessi0.reslt : ../src/xbessi0.f $(ILIB)
	-if [ -r ../bin/xbessi0.out ]; then \
	../bin/xbessi0.out > xbessi0.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessi0.out ../src/xbessi0.f $(ILIB) && \
	xbessi0.out > xbessi0.reslt && /bin/rm -f xbessi0.o xbessi0.out; \
	fi

xbessi1.reslt : ../src/xbessi1.f $(ILIB)
	-if [ -r ../bin/xbessi1.out ]; then \
	../bin/xbessi1.out > xbessi1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessi1.out ../src/xbessi1.f $(ILIB) && \
	xbessi1.out > xbessi1.reslt && /bin/rm -f xbessi1.o xbessi1.out; \
	fi

xbessik.reslt : ../src/xbessik.f $(ILIB)
	-if [ -r ../bin/xbessik.out ]; then \
	../bin/xbessik.out > xbessik.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessik.out ../src/xbessik.f $(ILIB) && \
	xbessik.out > xbessik.reslt && /bin/rm -f xbessik.o xbessik.out; \
	fi

xbessj.reslt : ../src/xbessj.f $(ILIB)
	-if [ -r ../bin/xbessj.out ]; then \
	../bin/xbessj.out > xbessj.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessj.out ../src/xbessj.f $(ILIB) && \
	xbessj.out > xbessj.reslt && /bin/rm -f xbessj.o xbessj.out; \
	fi

xbessj0.reslt : ../src/xbessj0.f $(ILIB)
	-if [ -r ../bin/xbessj0.out ]; then \
	../bin/xbessj0.out > xbessj0.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessj0.out ../src/xbessj0.f $(ILIB) && \
	xbessj0.out > xbessj0.reslt && /bin/rm -f xbessj0.o xbessj0.out; \
	fi

xbessj1.reslt : ../src/xbessj1.f $(ILIB)
	-if [ -r ../bin/xbessj1.out ]; then \
	../bin/xbessj1.out > xbessj1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessj1.out ../src/xbessj1.f $(ILIB) && \
	xbessj1.out > xbessj1.reslt && /bin/rm -f xbessj1.o xbessj1.out; \
	fi

xbessjy.reslt : ../src/xbessjy.f $(ILIB)
	-if [ -r ../bin/xbessjy.out ]; then \
	../bin/xbessjy.out > xbessjy.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessjy.out ../src/xbessjy.f $(ILIB) && \
	xbessjy.out > xbessjy.reslt && /bin/rm -f xbessjy.o xbessjy.out; \
	fi

xbessk.reslt : ../src/xbessk.f $(ILIB)
	-if [ -r ../bin/xbessk.out ]; then \
	../bin/xbessk.out > xbessk.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessk.out ../src/xbessk.f $(ILIB) && \
	xbessk.out > xbessk.reslt && /bin/rm -f xbessk.o xbessk.out; \
	fi

xbessk0.reslt : ../src/xbessk0.f $(ILIB)
	-if [ -r ../bin/xbessk0.out ]; then \
	../bin/xbessk0.out > xbessk0.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessk0.out ../src/xbessk0.f $(ILIB) && \
	xbessk0.out > xbessk0.reslt && /bin/rm -f xbessk0.o xbessk0.out; \
	fi

xbessk1.reslt : ../src/xbessk1.f $(ILIB)
	-if [ -r ../bin/xbessk1.out ]; then \
	../bin/xbessk1.out > xbessk1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessk1.out ../src/xbessk1.f $(ILIB) && \
	xbessk1.out > xbessk1.reslt && /bin/rm -f xbessk1.o xbessk1.out; \
	fi

xbessy.reslt : ../src/xbessy.f $(ILIB)
	-if [ -r ../bin/xbessy.out ]; then \
	../bin/xbessy.out > xbessy.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessy.out ../src/xbessy.f $(ILIB) && \
	xbessy.out > xbessy.reslt && /bin/rm -f xbessy.o xbessy.out; \
	fi

xbessy0.reslt : ../src/xbessy0.f $(ILIB)
	-if [ -r ../bin/xbessy0.out ]; then \
	../bin/xbessy0.out > xbessy0.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessy0.out ../src/xbessy0.f $(ILIB) && \
	xbessy0.out > xbessy0.reslt && /bin/rm -f xbessy0.o xbessy0.out; \
	fi

xbessy1.reslt : ../src/xbessy1.f $(ILIB)
	-if [ -r ../bin/xbessy1.out ]; then \
	../bin/xbessy1.out > xbessy1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbessy1.out ../src/xbessy1.f $(ILIB) && \
	xbessy1.out > xbessy1.reslt && /bin/rm -f xbessy1.o xbessy1.out; \
	fi

xbeta.reslt : ../src/xbeta.f $(ILIB)
	-if [ -r ../bin/xbeta.out ]; then \
	../bin/xbeta.out > xbeta.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbeta.out ../src/xbeta.f $(ILIB) && \
	xbeta.out > xbeta.reslt && /bin/rm -f xbeta.o xbeta.out; \
	fi

xbetai.reslt : ../src/xbetai.f $(ILIB)
	-if [ -r ../bin/xbetai.out ]; then \
	../bin/xbetai.out > xbetai.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbetai.out ../src/xbetai.f $(ILIB) && \
	xbetai.out > xbetai.reslt && /bin/rm -f xbetai.o xbetai.out; \
	fi

xbico.reslt : ../src/xbico.f $(ILIB)
	-if [ -r ../bin/xbico.out ]; then \
	../bin/xbico.out > xbico.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbico.out ../src/xbico.f $(ILIB) && \
	xbico.out > xbico.reslt && /bin/rm -f xbico.o xbico.out; \
	fi

xbnldev.reslt : ../src/xbnldev.f $(ILIB)
	-if [ -r ../bin/xbnldev.out ]; then \
	../bin/xbnldev.out < ../responses/xbnldev.resp > xbnldev.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbnldev.out ../src/xbnldev.f $(ILIB) && \
	xbnldev.out < ../responses/xbnldev.resp > xbnldev.reslt && \
	/bin/rm -f xbnldev.o xbnldev.out; \
	fi

xbrent.reslt : ../src/xbrent.f $(ILIB)
	-if [ -r ../bin/xbrent.out ]; then \
	../bin/xbrent.out > xbrent.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbrent.out ../src/xbrent.f $(ILIB) && \
	xbrent.out > xbrent.reslt && /bin/rm -f xbrent.o xbrent.out; \
	fi

xbroydn.reslt : ../src/xbroydn.f $(ILIB)
	-if [ -r ../bin/xbroydn.out ]; then \
	../bin/xbroydn.out > xbroydn.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbroydn.out ../src/xbroydn.f $(ILIB) && \
	xbroydn.out > xbroydn.reslt && /bin/rm -f xbroydn.o xbroydn.out; \
	fi

xbsstep.reslt : ../src/xbsstep.f $(ILIB)
	-if [ -r ../bin/xbsstep.out ]; then \
	../bin/xbsstep.out > xbsstep.reslt; \
	else \
	$(FC) $(FFLAGS) -o xbsstep.out ../src/xbsstep.f $(ILIB) && \
	xbsstep.out > xbsstep.reslt && /bin/rm -f xbsstep.o xbsstep.out; \
	fi

xcaldat.reslt : ../src/xcaldat.f $(ILIB)
	-if [ -r ../bin/xcaldat.out ]; then \
	../bin/xcaldat.out > xcaldat.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcaldat.out ../src/xcaldat.f $(ILIB) && \
	xcaldat.out > xcaldat.reslt && /bin/rm -f xcaldat.o xcaldat.out; \
	fi

xchder.reslt : ../src/xchder.f $(ILIB)
	-if [ -r ../bin/xchder.out ]; then \
	../bin/xchder.out < ../responses/xchder.resp > xchder.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchder.out ../src/xchder.f $(ILIB) && \
	xchder.out < ../responses/xchder.resp > xchder.reslt && \
	/bin/rm -f xchder.o xchder.out; \
	fi

xchebev.reslt : ../src/xchebev.f $(ILIB)
	-if [ -r ../bin/xchebev.out ]; then \
	../bin/xchebev.out < ../responses/xchebev.resp > xchebev.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchebev.out ../src/xchebev.f $(ILIB) && \
	xchebev.out < ../responses/xchebev.resp > xchebev.reslt && \
	/bin/rm -f xchebev.o xchebev.out; \
	fi

xchebft.reslt : ../src/xchebft.f $(ILIB)
	-if [ -r ../bin/xchebft.out ]; then \
	../bin/xchebft.out < ../responses/xchebft.resp > xchebft.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchebft.out ../src/xchebft.f $(ILIB) && \
	xchebft.out < ../responses/xchebft.resp > xchebft.reslt && \
	/bin/rm -f xchebft.o xchebft.out; \
	fi

xchebpc.reslt : ../src/xchebpc.f $(ILIB)
	-if [ -r ../bin/xchebpc.out ]; then \
	../bin/xchebpc.out < ../responses/xchebpc.resp > xchebpc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchebpc.out ../src/xchebpc.f $(ILIB) && \
	xchebpc.out < ../responses/xchebpc.resp > xchebpc.reslt && \
	/bin/rm -f xchebpc.o xchebpc.out; \
	fi

xchint.reslt : ../src/xchint.f $(ILIB)
	-if [ -r ../bin/xchint.out ]; then \
	../bin/xchint.out < ../responses/xchint.resp > xchint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchint.out ../src/xchint.f $(ILIB) && \
	xchint.out < ../responses/xchint.resp > xchint.reslt && \
	/bin/rm -f xchint.o xchint.out; \
	fi

xcholsl.reslt : ../src/xcholsl.f $(ILIB)
	-if [ -r ../bin/xcholsl.out ]; then \
	../bin/xcholsl.out > xcholsl.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcholsl.out ../src/xcholsl.f $(ILIB) && \
	xcholsl.out > xcholsl.reslt && /bin/rm -f xcholsl.o xcholsl.out; \
	fi

xchsone.reslt : ../src/xchsone.f $(ILIB)
	-if [ -r ../bin/xchsone.out ]; then \
	../bin/xchsone.out > xchsone.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchsone.out ../src/xchsone.f $(ILIB) && \
	xchsone.out > xchsone.reslt && /bin/rm -f xchsone.o xchsone.out; \
	fi

xchstwo.reslt : ../src/xchstwo.f $(ILIB)
	-if [ -r ../bin/xchstwo.out ]; then \
	../bin/xchstwo.out > xchstwo.reslt; \
	else \
	$(FC) $(FFLAGS) -o xchstwo.out ../src/xchstwo.f $(ILIB) && \
	xchstwo.out > xchstwo.reslt && /bin/rm -f xchstwo.o xchstwo.out; \
	fi

xcisi.reslt : ../src/xcisi.f $(ILIB)
	-if [ -r ../bin/xcisi.out ]; then \
	../bin/xcisi.out > xcisi.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcisi.out ../src/xcisi.f $(ILIB) && \
	xcisi.out > xcisi.reslt && /bin/rm -f xcisi.o xcisi.out; \
	fi

xcntab1.reslt : ../src/xcntab1.f $(ILIB)
	-if [ -r ../bin/xcntab1.out ]; then \
	../bin/xcntab1.out > xcntab1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcntab1.out ../src/xcntab1.f $(ILIB) && \
	xcntab1.out > xcntab1.reslt && /bin/rm -f xcntab1.o xcntab1.out; \
	fi

xcntab2.reslt : ../src/xcntab2.f $(ILIB)
	-if [ -r ../bin/xcntab2.out ]; then \
	../bin/xcntab2.out > xcntab2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcntab2.out ../src/xcntab2.f $(ILIB) && \
	xcntab2.out > xcntab2.reslt && /bin/rm -f xcntab2.o xcntab2.out; \
	fi

xconvlv.reslt : ../src/xconvlv.f $(ILIB)
	-if [ -r ../bin/xconvlv.out ]; then \
	../bin/xconvlv.out > xconvlv.reslt; \
	else \
	$(FC) $(FFLAGS) -o xconvlv.out ../src/xconvlv.f $(ILIB) && \
	xconvlv.out > xconvlv.reslt && /bin/rm -f xconvlv.o xconvlv.out; \
	fi

xcorrel.reslt : ../src/xcorrel.f $(ILIB)
	-if [ -r ../bin/xcorrel.out ]; then \
	../bin/xcorrel.out > xcorrel.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcorrel.out ../src/xcorrel.f $(ILIB) && \
	xcorrel.out > xcorrel.reslt && /bin/rm -f xcorrel.o xcorrel.out; \
	fi

xcosft1.reslt : ../src/xcosft1.f $(ILIB)
	-if [ -r ../bin/xcosft1.out ]; then \
	../bin/xcosft1.out < ../responses/xcosft1.resp > xcosft1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcosft1.out ../src/xcosft1.f $(ILIB) && \
	xcosft1.out < ../responses/xcosft1.resp > xcosft1.reslt && \
	/bin/rm -f xcosft1.o xcosft1.out; \
	fi

xcosft2.reslt : ../src/xcosft2.f $(ILIB)
	-if [ -r ../bin/xcosft2.out ]; then \
	../bin/xcosft2.out < ../responses/xcosft2.resp > xcosft2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcosft2.out ../src/xcosft2.f $(ILIB) && \
	xcosft2.out < ../responses/xcosft2.resp > xcosft2.reslt && \
	/bin/rm -f xcosft2.o xcosft2.out; \
	fi

xcovsrt.reslt : ../src/xcovsrt.f $(ILIB)
	-if [ -r ../bin/xcovsrt.out ]; then \
	../bin/xcovsrt.out < ../responses/xcovsrt.resp > xcovsrt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcovsrt.out ../src/xcovsrt.f $(ILIB) && \
	xcovsrt.out < ../responses/xcovsrt.resp > xcovsrt.reslt && \
	/bin/rm -f xcovsrt.o xcovsrt.out; \
	fi

xcrank.reslt : ../src/xcrank.f $(ILIB)
	-if [ -r ../bin/xcrank.out ]; then \
	../bin/xcrank.out < ../responses/xcrank.resp > xcrank.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcrank.out ../src/xcrank.f $(ILIB) && \
	xcrank.out < ../responses/xcrank.resp > xcrank.reslt && \
	/bin/rm -f xcrank.o xcrank.out; \
	fi

xcyclic.reslt : ../src/xcyclic.f $(ILIB)
	-if [ -r ../bin/xcyclic.out ]; then \
	../bin/xcyclic.out > xcyclic.reslt; \
	else \
	$(FC) $(FFLAGS) -o xcyclic.out ../src/xcyclic.f $(ILIB) && \
	xcyclic.out > xcyclic.reslt && /bin/rm -f xcyclic.o xcyclic.out; \
	fi

xdawson.reslt : ../src/xdawson.f $(ILIB)
	-if [ -r ../bin/xdawson.out ]; then \
	../bin/xdawson.out > xdawson.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdawson.out ../src/xdawson.f $(ILIB) && \
	xdawson.out > xdawson.reslt && /bin/rm -f xdawson.o xdawson.out; \
	fi

xdbrent.reslt : ../src/xdbrent.f $(ILIB)
	-if [ -r ../bin/xdbrent.out ]; then \
	../bin/xdbrent.out > xdbrent.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdbrent.out ../src/xdbrent.f $(ILIB) && \
	xdbrent.out > xdbrent.reslt && /bin/rm -f xdbrent.o xdbrent.out; \
	fi

xddpoly.reslt : ../src/xddpoly.f $(ILIB)
	-if [ -r ../bin/xddpoly.out ]; then \
	../bin/xddpoly.out < ../responses/xddpoly.resp > xddpoly.reslt; \
	else \
	$(FC) $(FFLAGS) -o xddpoly.out ../src/xddpoly.f $(ILIB) && \
	xddpoly.out < ../responses/xddpoly.resp > xddpoly.reslt && \
	/bin/rm -f xddpoly.o xddpoly.out; \
	fi

xdecchk.reslt : ../src/xdecchk.f $(ILIB)
	-if [ -r ../bin/xdecchk.out ]; then \
	../bin/xdecchk.out < ../responses/xdecchk.resp > xdecchk.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdecchk.out ../src/xdecchk.f $(ILIB) && \
	xdecchk.out < ../responses/xdecchk.resp > xdecchk.reslt && \
	/bin/rm -f xdecchk.o xdecchk.out; \
	fi

xdf1dim.reslt : ../src/xdf1dim.f $(ILIB)
	-if [ -r ../bin/xdf1dim.out ]; then \
	../bin/xdf1dim.out < ../responses/xdf1dim.resp > xdf1dim.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdf1dim.out ../src/xdf1dim.f $(ILIB) && \
	xdf1dim.out < ../responses/xdf1dim.resp > xdf1dim.reslt && \
	/bin/rm -f xdf1dim.o xdf1dim.out; \
	fi

xdfpmin.reslt : ../src/xdfpmin.f $(ILIB)
	-if [ -r ../bin/xdfpmin.out ]; then \
	../bin/xdfpmin.out > xdfpmin.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdfpmin.out ../src/xdfpmin.f $(ILIB) && \
	xdfpmin.out > xdfpmin.reslt && /bin/rm -f xdfpmin.o xdfpmin.out; \
	fi

xdfridr.reslt : ../src/xdfridr.f $(ILIB)
	-if [ -r ../bin/xdfridr.out ]; then \
	../bin/xdfridr.out < ../responses/xdfridr.resp > xdfridr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdfridr.out ../src/xdfridr.f $(ILIB) && \
	xdfridr.out < ../responses/xdfridr.resp > xdfridr.reslt && \
	/bin/rm -f xdfridr.o xdfridr.out; \
	fi

xdftint.reslt : ../src/xdftint.f $(ILIB)
	-if [ -r ../bin/xdftint.out ]; then \
	../bin/xdftint.out < ../responses/xdftint.resp > xdftint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xdftint.out ../src/xdftint.f $(ILIB) && \
	xdftint.out < ../responses/xdftint.resp > xdftint.reslt && \
	/bin/rm -f xdftint.o xdftint.out; \
	fi

xeclass.reslt : ../src/xeclass.f $(ILIB)
	-if [ -r ../bin/xeclass.out ]; then \
	../bin/xeclass.out > xeclass.reslt; \
	else \
	$(FC) $(FFLAGS) -o xeclass.out ../src/xeclass.f $(ILIB) && \
	xeclass.out > xeclass.reslt && /bin/rm -f xeclass.o xeclass.out; \
	fi

xeclazz.reslt : ../src/xeclazz.f $(ILIB)
	-if [ -r ../bin/xeclazz.out ]; then \
	../bin/xeclazz.out > xeclazz.reslt; \
	else \
	$(FC) $(FFLAGS) -o xeclazz.out ../src/xeclazz.f $(ILIB) && \
	xeclazz.out > xeclazz.reslt && /bin/rm -f xeclazz.o xeclazz.out; \
	fi

xei.reslt : ../src/xei.f $(ILIB)
	-if [ -r ../bin/xei.out ]; then \
	../bin/xei.out > xei.reslt; \
	else \
	$(FC) $(FFLAGS) -o xei.out ../src/xei.f $(ILIB) && \
	xei.out > xei.reslt && /bin/rm -f xei.o xei.out; \
	fi

xeigsrt.reslt : ../src/xeigsrt.f $(ILIB)
	-if [ -r ../bin/xeigsrt.out ]; then \
	../bin/xeigsrt.out > xeigsrt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xeigsrt.out ../src/xeigsrt.f $(ILIB) && \
	xeigsrt.out > xeigsrt.reslt && /bin/rm -f xeigsrt.o xeigsrt.out; \
	fi

xelle.reslt : ../src/xelle.f $(ILIB)
	-if [ -r ../bin/xelle.out ]; then \
	../bin/xelle.out > xelle.reslt; \
	else \
	$(FC) $(FFLAGS) -o xelle.out ../src/xelle.f $(ILIB) && \
	xelle.out > xelle.reslt && /bin/rm -f xelle.o xelle.out; \
	fi

xellf.reslt : ../src/xellf.f $(ILIB)
	-if [ -r ../bin/xellf.out ]; then \
	../bin/xellf.out > xellf.reslt; \
	else \
	$(FC) $(FFLAGS) -o xellf.out ../src/xellf.f $(ILIB) && \
	xellf.out > xellf.reslt && /bin/rm -f xellf.o xellf.out; \
	fi

xellpi.reslt : ../src/xellpi.f $(ILIB)
	-if [ -r ../bin/xellpi.out ]; then \
	../bin/xellpi.out > xellpi.reslt; \
	else \
	$(FC) $(FFLAGS) -o xellpi.out ../src/xellpi.f $(ILIB) && \
	xellpi.out > xellpi.reslt && /bin/rm -f xellpi.o xellpi.out; \
	fi

xelmhes.reslt : ../src/xelmhes.f $(ILIB)
	-if [ -r ../bin/xelmhes.out ]; then \
	../bin/xelmhes.out > xelmhes.reslt; \
	else \
	$(FC) $(FFLAGS) -o xelmhes.out ../src/xelmhes.f $(ILIB) && \
	xelmhes.out > xelmhes.reslt && /bin/rm -f xelmhes.o xelmhes.out; \
	fi

xerf.reslt : ../src/xerf.f $(ILIB)
	-if [ -r ../bin/xerf.out ]; then \
	../bin/xerf.out > xerf.reslt; \
	else \
	$(FC) $(FFLAGS) -o xerf.out ../src/xerf.f $(ILIB) && \
	xerf.out > xerf.reslt && /bin/rm -f xerf.o xerf.out; \
	fi

xerfc.reslt : ../src/xerfc.f $(ILIB)
	-if [ -r ../bin/xerfc.out ]; then \
	../bin/xerfc.out > xerfc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xerfc.out ../src/xerfc.f $(ILIB) && \
	xerfc.out > xerfc.reslt && /bin/rm -f xerfc.o xerfc.out; \
	fi

xerfcc.reslt : ../src/xerfcc.f $(ILIB)
	-if [ -r ../bin/xerfcc.out ]; then \
	../bin/xerfcc.out > xerfcc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xerfcc.out ../src/xerfcc.f $(ILIB) && \
	xerfcc.out > xerfcc.reslt && /bin/rm -f xerfcc.o xerfcc.out; \
	fi

xeulsum.reslt : ../src/xeulsum.f $(ILIB)
	-if [ -r ../bin/xeulsum.out ]; then \
	../bin/xeulsum.out < ../responses/xeulsum.resp > xeulsum.reslt; \
	else \
	$(FC) $(FFLAGS) -o xeulsum.out ../src/xeulsum.f $(ILIB) && \
	xeulsum.out < ../responses/xeulsum.resp > xeulsum.reslt && \
	/bin/rm -f xeulsum.o xeulsum.out; \
	fi

xevlmem.reslt : ../src/xevlmem.f $(ILIB)
	-if [ -r ../bin/xevlmem.out ]; then \
	../bin/xevlmem.out > xevlmem.reslt; \
	else \
	$(FC) $(FFLAGS) -o xevlmem.out ../src/xevlmem.f $(ILIB) && \
	xevlmem.out > xevlmem.reslt && /bin/rm -f xevlmem.o xevlmem.out; \
	fi

xexpdev.reslt : ../src/xexpdev.f $(ILIB)
	-if [ -r ../bin/xexpdev.out ]; then \
	../bin/xexpdev.out > xexpdev.reslt; \
	else \
	$(FC) $(FFLAGS) -o xexpdev.out ../src/xexpdev.f $(ILIB) && \
	xexpdev.out > xexpdev.reslt && /bin/rm -f xexpdev.o xexpdev.out; \
	fi

xexpint.reslt : ../src/xexpint.f $(ILIB)
	-if [ -r ../bin/xexpint.out ]; then \
	../bin/xexpint.out > xexpint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xexpint.out ../src/xexpint.f $(ILIB) && \
	xexpint.out > xexpint.reslt && /bin/rm -f xexpint.o xexpint.out; \
	fi

xf1dim.reslt : ../src/xf1dim.f $(ILIB)
	-if [ -r ../bin/xf1dim.out ]; then \
	../bin/xf1dim.out < ../responses/xf1dim.resp > xf1dim.reslt; \
	else \
	$(FC) $(FFLAGS) -o xf1dim.out ../src/xf1dim.f $(ILIB) && \
	xf1dim.out < ../responses/xf1dim.resp > xf1dim.reslt && \
	/bin/rm -f xf1dim.o xf1dim.out; \
	fi

xfactln.reslt : ../src/xfactln.f $(ILIB)
	-if [ -r ../bin/xfactln.out ]; then \
	../bin/xfactln.out > xfactln.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfactln.out ../src/xfactln.f $(ILIB) && \
	xfactln.out > xfactln.reslt && /bin/rm -f xfactln.o xfactln.out; \
	fi

xfactrl.reslt : ../src/xfactrl.f $(ILIB)
	-if [ -r ../bin/xfactrl.out ]; then \
	../bin/xfactrl.out > xfactrl.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfactrl.out ../src/xfactrl.f $(ILIB) && \
	xfactrl.out > xfactrl.reslt && /bin/rm -f xfactrl.o xfactrl.out; \
	fi

xfasper.reslt : ../src/xfasper.f $(ILIB)
	-if [ -r ../bin/xfasper.out ]; then \
	../bin/xfasper.out > xfasper.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfasper.out ../src/xfasper.f $(ILIB) && \
	xfasper.out > xfasper.reslt && /bin/rm -f xfasper.o xfasper.out; \
	fi

xfgauss.reslt : ../src/xfgauss.f $(ILIB)
	-if [ -r ../bin/xfgauss.out ]; then \
	../bin/xfgauss.out > xfgauss.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfgauss.out ../src/xfgauss.f $(ILIB) && \
	xfgauss.out > xfgauss.reslt && /bin/rm -f xfgauss.o xfgauss.out; \
	fi

xfit.reslt : ../src/xfit.f $(ILIB)
	-if [ -r ../bin/xfit.out ]; then \
	../bin/xfit.out > xfit.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfit.out ../src/xfit.f $(ILIB) && \
	xfit.out > xfit.reslt && /bin/rm -f xfit.o xfit.out; \
	fi

xfitexy.reslt : ../src/xfitexy.f $(ILIB)
	-if [ -r ../bin/xfitexy.out ]; then \
	../bin/xfitexy.out > xfitexy.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfitexy.out ../src/xfitexy.f $(ILIB) && \
	xfitexy.out > xfitexy.reslt && /bin/rm -f xfitexy.o xfitexy.out; \
	fi

xfixrts.reslt : ../src/xfixrts.f $(ILIB)
	-if [ -r ../bin/xfixrts.out ]; then \
	../bin/xfixrts.out > xfixrts.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfixrts.out ../src/xfixrts.f $(ILIB) && \
	xfixrts.out > xfixrts.reslt && /bin/rm -f xfixrts.o xfixrts.out; \
	fi

xfleg.reslt : ../src/xfleg.f $(ILIB)
	-if [ -r ../bin/xfleg.out ]; then \
	../bin/xfleg.out > xfleg.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfleg.out ../src/xfleg.f $(ILIB) && \
	xfleg.out > xfleg.reslt && /bin/rm -f xfleg.o xfleg.out; \
	fi

xflmoon.reslt : ../src/xflmoon.f $(ILIB)
	-if [ -r ../bin/xflmoon.out ]; then \
	../bin/xflmoon.out < ../responses/xflmoon.resp > xflmoon.reslt; \
	else \
	$(FC) $(FFLAGS) -o xflmoon.out ../src/xflmoon.f $(ILIB) && \
	xflmoon.out < ../responses/xflmoon.resp > xflmoon.reslt && \
	/bin/rm -f xflmoon.o xflmoon.out; \
	fi

xfour1.reslt : ../src/xfour1.f $(ILIB)
	-if [ -r ../bin/xfour1.out ]; then \
	../bin/xfour1.out < ../responses/xfour1.resp > xfour1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfour1.out ../src/xfour1.f $(ILIB) && \
	xfour1.out < ../responses/xfour1.resp > xfour1.reslt && \
	/bin/rm -f xfour1.o xfour1.out; \
	fi

xfourfs.reslt : ../src/xfourfs.f $(ILIB)
	-if [ -r ../bin/xfourfs.out ]; then \
	../bin/xfourfs.out > xfourfs.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfourfs.out ../src/xfourfs.f $(ILIB) && \
	xfourfs.out > xfourfs.reslt && /bin/rm -f xfourfs.o xfourfs.out; \
	fi

xfourn.reslt : ../src/xfourn.f $(ILIB)
	-if [ -r ../bin/xfourn.out ]; then \
	../bin/xfourn.out > xfourn.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfourn.out ../src/xfourn.f $(ILIB) && \
	xfourn.out > xfourn.reslt && /bin/rm -f xfourn.o xfourn.out; \
	fi

xfpoly.reslt : ../src/xfpoly.f $(ILIB)
	-if [ -r ../bin/xfpoly.out ]; then \
	../bin/xfpoly.out > xfpoly.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfpoly.out ../src/xfpoly.f $(ILIB) && \
	xfpoly.out > xfpoly.reslt && /bin/rm -f xfpoly.o xfpoly.out; \
	fi

xfred2.reslt : ../src/xfred2.f $(ILIB)
	-if [ -r ../bin/xfred2.out ]; then \
	../bin/xfred2.out > xfred2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfred2.out ../src/xfred2.f $(ILIB) && \
	xfred2.out > xfred2.reslt && /bin/rm -f xfred2.o xfred2.out; \
	fi

xfredin.reslt : ../src/xfredin.f $(ILIB)
	-if [ -r ../bin/xfredin.out ]; then \
	../bin/xfredin.out < ../responses/xfredin.resp > xfredin.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfredin.out ../src/xfredin.f $(ILIB) && \
	xfredin.out < ../responses/xfredin.resp > xfredin.reslt && \
	/bin/rm -f xfredin.o xfredin.out; \
	fi

xfrenel.reslt : ../src/xfrenel.f $(ILIB)
	-if [ -r ../bin/xfrenel.out ]; then \
	../bin/xfrenel.out > xfrenel.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfrenel.out ../src/xfrenel.f $(ILIB) && \
	xfrenel.out > xfrenel.reslt && /bin/rm -f xfrenel.o xfrenel.out; \
	fi

xfrprmn.reslt : ../src/xfrprmn.f $(ILIB)
	-if [ -r ../bin/xfrprmn.out ]; then \
	../bin/xfrprmn.out > xfrprmn.reslt; \
	else \
	$(FC) $(FFLAGS) -o xfrprmn.out ../src/xfrprmn.f $(ILIB) && \
	xfrprmn.out > xfrprmn.reslt && /bin/rm -f xfrprmn.o xfrprmn.out; \
	fi

xftest.reslt : ../src/xftest.f $(ILIB)
	-if [ -r ../bin/xftest.out ]; then \
	../bin/xftest.out > xftest.reslt; \
	else \
	$(FC) $(FFLAGS) -o xftest.out ../src/xftest.f $(ILIB) && \
	xftest.out > xftest.reslt && /bin/rm -f xftest.o xftest.out; \
	fi

xgamdev.reslt : ../src/xgamdev.f $(ILIB)
	-if [ -r ../bin/xgamdev.out ]; then \
	../bin/xgamdev.out < ../responses/xgamdev.resp > xgamdev.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgamdev.out ../src/xgamdev.f $(ILIB) && \
	xgamdev.out < ../responses/xgamdev.resp > xgamdev.reslt && \
	/bin/rm -f xgamdev.o xgamdev.out; \
	fi

xgammln.reslt : ../src/xgammln.f $(ILIB)
	-if [ -r ../bin/xgammln.out ]; then \
	../bin/xgammln.out > xgammln.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgammln.out ../src/xgammln.f $(ILIB) && \
	xgammln.out > xgammln.reslt && /bin/rm -f xgammln.o xgammln.out; \
	fi

xgammp.reslt : ../src/xgammp.f $(ILIB)
	-if [ -r ../bin/xgammp.out ]; then \
	../bin/xgammp.out > xgammp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgammp.out ../src/xgammp.f $(ILIB) && \
	xgammp.out > xgammp.reslt && /bin/rm -f xgammp.o xgammp.out; \
	fi

xgammq.reslt : ../src/xgammq.f $(ILIB)
	-if [ -r ../bin/xgammq.out ]; then \
	../bin/xgammq.out > xgammq.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgammq.out ../src/xgammq.f $(ILIB) && \
	xgammq.out > xgammq.reslt && /bin/rm -f xgammq.o xgammq.out; \
	fi

xgasdev.reslt : ../src/xgasdev.f $(ILIB)
	-if [ -r ../bin/xgasdev.out ]; then \
	../bin/xgasdev.out > xgasdev.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgasdev.out ../src/xgasdev.f $(ILIB) && \
	xgasdev.out > xgasdev.reslt && /bin/rm -f xgasdev.o xgasdev.out; \
	fi

xgaucof.reslt : ../src/xgaucof.f $(ILIB)
	-if [ -r ../bin/xgaucof.out ]; then \
	../bin/xgaucof.out < ../responses/xgaucof.resp > xgaucof.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgaucof.out ../src/xgaucof.f $(ILIB) && \
	xgaucof.out < ../responses/xgaucof.resp > xgaucof.reslt && \
	/bin/rm -f xgaucof.o xgaucof.out; \
	fi

xgauher.reslt : ../src/xgauher.f $(ILIB)
	-if [ -r ../bin/xgauher.out ]; then \
	../bin/xgauher.out < ../responses/xgauher.resp > xgauher.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgauher.out ../src/xgauher.f $(ILIB) && \
	xgauher.out < ../responses/xgauher.resp > xgauher.reslt && \
	/bin/rm -f xgauher.o xgauher.out; \
	fi

xgaujac.reslt : ../src/xgaujac.f $(ILIB)
	-if [ -r ../bin/xgaujac.out ]; then \
	../bin/xgaujac.out < ../responses/xgaujac.resp > xgaujac.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgaujac.out ../src/xgaujac.f $(ILIB) && \
	xgaujac.out < ../responses/xgaujac.resp > xgaujac.reslt && \
	/bin/rm -f xgaujac.o xgaujac.out; \
	fi

xgaulag.reslt : ../src/xgaulag.f $(ILIB)
	-if [ -r ../bin/xgaulag.out ]; then \
	../bin/xgaulag.out < ../responses/xgaulag.resp > xgaulag.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgaulag.out ../src/xgaulag.f $(ILIB) && \
	xgaulag.out < ../responses/xgaulag.resp > xgaulag.reslt && \
	/bin/rm -f xgaulag.o xgaulag.out; \
	fi

xgauleg.reslt : ../src/xgauleg.f $(ILIB)
	-if [ -r ../bin/xgauleg.out ]; then \
	../bin/xgauleg.out > xgauleg.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgauleg.out ../src/xgauleg.f $(ILIB) && \
	xgauleg.out > xgauleg.reslt && /bin/rm -f xgauleg.o xgauleg.out; \
	fi

xgaussj.reslt : ../src/xgaussj.f $(ILIB)
	-if [ -r ../bin/xgaussj.out ]; then \
	../bin/xgaussj.out < ../responses/xgaussj.resp > xgaussj.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgaussj.out ../src/xgaussj.f $(ILIB) && \
	xgaussj.out < ../responses/xgaussj.resp > xgaussj.reslt && \
	/bin/rm -f xgaussj.o xgaussj.out; \
	fi

xgcf.reslt : ../src/xgcf.f $(ILIB)
	-if [ -r ../bin/xgcf.out ]; then \
	../bin/xgcf.out > xgcf.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgcf.out ../src/xgcf.f $(ILIB) && \
	xgcf.out > xgcf.reslt && /bin/rm -f xgcf.o xgcf.out; \
	fi

xgolden.reslt : ../src/xgolden.f $(ILIB)
	-if [ -r ../bin/xgolden.out ]; then \
	../bin/xgolden.out > xgolden.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgolden.out ../src/xgolden.f $(ILIB) && \
	xgolden.out > xgolden.reslt && /bin/rm -f xgolden.o xgolden.out; \
	fi

xgser.reslt : ../src/xgser.f $(ILIB)
	-if [ -r ../bin/xgser.out ]; then \
	../bin/xgser.out > xgser.reslt; \
	else \
	$(FC) $(FFLAGS) -o xgser.out ../src/xgser.f $(ILIB) && \
	xgser.out > xgser.reslt && /bin/rm -f xgser.o xgser.out; \
	fi

xhexout.reslt : ../src/xhexout.f $(ILIB)
	-if [ -r ../bin/xhexout.out ]; then \
	../bin/xhexout.out > xhexout.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhexout.out ../src/xhexout.f $(ILIB) && \
	xhexout.out > xhexout.reslt && /bin/rm -f xhexout.o xhexout.out; \
	fi

xhpsel.reslt : ../src/xhpsel.f $(ILIB)
	-if [ -r ../bin/xhpsel.out ]; then \
	../bin/xhpsel.out < ../responses/xhpsel.resp > xhpsel.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhpsel.out ../src/xhpsel.f $(ILIB) && \
	xhpsel.out < ../responses/xhpsel.resp > xhpsel.reslt && \
	/bin/rm -f xhpsel.o xhpsel.out; \
	fi

xhpsort.reslt : ../src/xhpsort.f $(ILIB)
	-if [ -r ../bin/xhpsort.out ]; then \
	../bin/xhpsort.out > xhpsort.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhpsort.out ../src/xhpsort.f $(ILIB) && \
	xhpsort.out > xhpsort.reslt && /bin/rm -f xhpsort.o xhpsort.out; \
	fi

xhqr.reslt : ../src/xhqr.f $(ILIB)
	-if [ -r ../bin/xhqr.out ]; then \
	../bin/xhqr.out > xhqr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhqr.out ../src/xhqr.f $(ILIB) && \
	xhqr.out > xhqr.reslt && /bin/rm -f xhqr.o xhqr.out; \
	fi

xhuffman.reslt : ../src/xhuffman.f $(ILIB)
	-if [ -r ../bin/xhuffman.out ]; then \
	../bin/xhuffman.out < ../responses/xhuffman.resp > xhuffman.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhuffman.out ../src/xhuffman.f $(ILIB) && \
	xhuffman.out < ../responses/xhuffman.resp > xhuffman.reslt && \
	/bin/rm -f xhuffman.o xhuffman.out; \
	fi

xhunt.reslt : ../src/xhunt.f $(ILIB)
	-if [ -r ../bin/xhunt.out ]; then \
	../bin/xhunt.out > xhunt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhunt.out ../src/xhunt.f $(ILIB) && \
	xhunt.out > xhunt.reslt && /bin/rm -f xhunt.o xhunt.out; \
	fi

xhypgeo.reslt : ../src/xhypgeo.f $(ILIB)
	-if [ -r ../bin/xhypgeo.out ]; then \
	../bin/xhypgeo.out < ../responses/xhypgeo.resp > xhypgeo.reslt; \
	else \
	$(FC) $(FFLAGS) -o xhypgeo.out ../src/xhypgeo.f $(ILIB) && \
	xhypgeo.out < ../responses/xhypgeo.resp > xhypgeo.reslt && \
	/bin/rm -f xhypgeo.o xhypgeo.out; \
	fi

xicrc.reslt : ../src/xicrc.f $(ILIB)
	-if [ -r ../bin/xicrc.out ]; then \
	../bin/xicrc.out < ../responses/xicrc.resp > xicrc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xicrc.out ../src/xicrc.f $(ILIB) && \
	xicrc.out < ../responses/xicrc.resp > xicrc.reslt && \
	/bin/rm -f xicrc.o xicrc.out; \
	fi

xigray.reslt : ../src/xigray.f $(ILIB)
	-if [ -r ../bin/xigray.out ]; then \
	../bin/xigray.out < ../responses/xigray.resp > xigray.reslt; \
	else \
	$(FC) $(FFLAGS) -o xigray.out ../src/xigray.f $(ILIB) && \
	xigray.out < ../responses/xigray.resp > xigray.reslt && \
	/bin/rm -f xigray.o xigray.out; \
	fi

xindexx.reslt : ../src/xindexx.f $(ILIB)
	-if [ -r ../bin/xindexx.out ]; then \
	../bin/xindexx.out > xindexx.reslt; \
	else \
	$(FC) $(FFLAGS) -o xindexx.out ../src/xindexx.f $(ILIB) && \
	xindexx.out > xindexx.reslt && /bin/rm -f xindexx.o xindexx.out; \
	fi

xirbit1.reslt : ../src/xirbit1.f $(ILIB)
	-if [ -r ../bin/xirbit1.out ]; then \
	../bin/xirbit1.out > xirbit1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xirbit1.out ../src/xirbit1.f $(ILIB) && \
	xirbit1.out > xirbit1.reslt && /bin/rm -f xirbit1.o xirbit1.out; \
	fi

xirbit2.reslt : ../src/xirbit2.f $(ILIB)
	-if [ -r ../bin/xirbit2.out ]; then \
	../bin/xirbit2.out > xirbit2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xirbit2.out ../src/xirbit2.f $(ILIB) && \
	xirbit2.out > xirbit2.reslt && /bin/rm -f xirbit2.o xirbit2.out; \
	fi

xjacobi.reslt : ../src/xjacobi.f $(ILIB)
	-if [ -r ../bin/xjacobi.out ]; then \
	../bin/xjacobi.out < ../responses/xjacobi.resp > xjacobi.reslt; \
	else \
	$(FC) $(FFLAGS) -o xjacobi.out ../src/xjacobi.f $(ILIB) && \
	xjacobi.out < ../responses/xjacobi.resp > xjacobi.reslt && \
	/bin/rm -f xjacobi.o xjacobi.out; \
	fi

xjulday.reslt : ../src/xjulday.f $(ILIB)
	-if [ -r ../bin/xjulday.out ]; then \
	../bin/xjulday.out < ../responses/xjulday.resp > xjulday.reslt; \
	else \
	$(FC) $(FFLAGS) -o xjulday.out ../src/xjulday.f $(ILIB) && \
	xjulday.out < ../responses/xjulday.resp > xjulday.reslt && \
	/bin/rm -f xjulday.o xjulday.out; \
	fi

xkendl1.reslt : ../src/xkendl1.f $(ILIB)
	-if [ -r ../bin/xkendl1.out ]; then \
	../bin/xkendl1.out > xkendl1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xkendl1.out ../src/xkendl1.f $(ILIB) && \
	xkendl1.out > xkendl1.reslt && /bin/rm -f xkendl1.o xkendl1.out; \
	fi

xkendl2.reslt : ../src/xkendl2.f $(ILIB)
	-if [ -r ../bin/xkendl2.out ]; then \
	../bin/xkendl2.out < ../responses/xkendl2.resp > xkendl2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xkendl2.out ../src/xkendl2.f $(ILIB) && \
	xkendl2.out < ../responses/xkendl2.resp > xkendl2.reslt && \
	/bin/rm -f xkendl2.o xkendl2.out; \
	fi

xks2d1s.reslt : ../src/xks2d1s.f $(ILIB)
	-if [ -r ../bin/xks2d1s.out ]; then \
	../bin/xks2d1s.out < ../responses/xks2d1s.resp > xks2d1s.reslt; \
	else \
	$(FC) $(FFLAGS) -o xks2d1s.out ../src/xks2d1s.f $(ILIB) && \
	xks2d1s.out < ../responses/xks2d1s.resp > xks2d1s.reslt && \
	/bin/rm -f xks2d1s.o xks2d1s.out; \
	fi

xks2d2s.reslt : ../src/xks2d2s.f $(ILIB)
	-if [ -r ../bin/xks2d2s.out ]; then \
	../bin/xks2d2s.out < ../responses/xks2d2s.resp > xks2d2s.reslt; \
	else \
	$(FC) $(FFLAGS) -o xks2d2s.out ../src/xks2d2s.f $(ILIB) && \
	xks2d2s.out < ../responses/xks2d2s.resp > xks2d2s.reslt && \
	/bin/rm -f xks2d2s.o xks2d2s.out; \
	fi

xksone.reslt : ../src/xksone.f $(ILIB)
	-if [ -r ../bin/xksone.out ]; then \
	../bin/xksone.out > xksone.reslt; \
	else \
	$(FC) $(FFLAGS) -o xksone.out ../src/xksone.f $(ILIB) && \
	xksone.out > xksone.reslt && /bin/rm -f xksone.o xksone.out; \
	fi

xkstwo.reslt : ../src/xkstwo.f $(ILIB)
	-if [ -r ../bin/xkstwo.out ]; then \
	../bin/xkstwo.out > xkstwo.reslt; \
	else \
	$(FC) $(FFLAGS) -o xkstwo.out ../src/xkstwo.f $(ILIB) && \
	xkstwo.out > xkstwo.reslt && /bin/rm -f xkstwo.o xkstwo.out; \
	fi

xlaguer.reslt : ../src/xlaguer.f $(ILIB)
	-if [ -r ../bin/xlaguer.out ]; then \
	../bin/xlaguer.out > xlaguer.reslt; \
	else \
	$(FC) $(FFLAGS) -o xlaguer.out ../src/xlaguer.f $(ILIB) && \
	xlaguer.out > xlaguer.reslt && /bin/rm -f xlaguer.o xlaguer.out; \
	fi

xlfit.reslt : ../src/xlfit.f $(ILIB)
	-if [ -r ../bin/xlfit.out ]; then \
	../bin/xlfit.out < ../responses/xlfit.resp > xlfit.reslt; \
	else \
	$(FC) $(FFLAGS) -o xlfit.out ../src/xlfit.f $(ILIB) && \
	xlfit.out < ../responses/xlfit.resp > xlfit.reslt && \
	/bin/rm -f xlfit.o xlfit.out; \
	fi

xlinbcg.reslt : ../src/xlinbcg.f $(ILIB)
	-if [ -r ../bin/xlinbcg.out ]; then \
	../bin/xlinbcg.out < ../responses/xlinbcg.resp > xlinbcg.reslt; \
	else \
	$(FC) $(FFLAGS) -o xlinbcg.out ../src/xlinbcg.f $(ILIB) && \
	xlinbcg.out < ../responses/xlinbcg.resp > xlinbcg.reslt && \
	/bin/rm -f xlinbcg.o xlinbcg.out; \
	fi

xlinmin.reslt : ../src/xlinmin.f $(ILIB)
	-if [ -r ../bin/xlinmin.out ]; then \
	../bin/xlinmin.out > xlinmin.reslt; \
	else \
	$(FC) $(FFLAGS) -o xlinmin.out ../src/xlinmin.f $(ILIB) && \
	xlinmin.out > xlinmin.reslt && /bin/rm -f xlinmin.o xlinmin.out; \
	fi

xlocate.reslt : ../src/xlocate.f $(ILIB)
	-if [ -r ../bin/xlocate.out ]; then \
	../bin/xlocate.out > xlocate.reslt; \
	else \
	$(FC) $(FFLAGS) -o xlocate.out ../src/xlocate.f $(ILIB) && \
	xlocate.out > xlocate.reslt && /bin/rm -f xlocate.o xlocate.out; \
	fi

xlubksb.reslt : ../src/xlubksb.f $(ILIB)
	-if [ -r ../bin/xlubksb.out ]; then \
	../bin/xlubksb.out < ../responses/xlubksb.resp > xlubksb.reslt; \
	else \
	$(FC) $(FFLAGS) -o xlubksb.out ../src/xlubksb.f $(ILIB) && \
	xlubksb.out < ../responses/xlubksb.resp > xlubksb.reslt && \
	/bin/rm -f xlubksb.o xlubksb.out; \
	fi

xludcmp.reslt : ../src/xludcmp.f $(ILIB)
	-if [ -r ../bin/xludcmp.out ]; then \
	../bin/xludcmp.out < ../responses/xludcmp.resp > xludcmp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xludcmp.out ../src/xludcmp.f $(ILIB) && \
	xludcmp.out < ../responses/xludcmp.resp > xludcmp.reslt && \
	/bin/rm -f xludcmp.o xludcmp.out; \
	fi

xmachar.reslt : ../src/xmachar.f $(ILIB)
	-if [ -r ../bin/xmachar.out ]; then \
	../bin/xmachar.out > xmachar.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmachar.out ../src/xmachar.f $(ILIB) && \
	xmachar.out > xmachar.reslt && /bin/rm -f xmachar.o xmachar.out; \
	fi

xmedfit.reslt : ../src/xmedfit.f $(ILIB)
	-if [ -r ../bin/xmedfit.out ]; then \
	../bin/xmedfit.out > xmedfit.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmedfit.out ../src/xmedfit.f $(ILIB) && \
	xmedfit.out > xmedfit.reslt && /bin/rm -f xmedfit.o xmedfit.out; \
	fi

xmemcof.reslt : ../src/xmemcof.f $(ILIB)
	-if [ -r ../bin/xmemcof.out ]; then \
	../bin/xmemcof.out > xmemcof.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmemcof.out ../src/xmemcof.f $(ILIB) && \
	xmemcof.out > xmemcof.reslt && /bin/rm -f xmemcof.o xmemcof.out; \
	fi

xmgfas.reslt : ../src/xmgfas.f $(ILIB)
	-if [ -r ../bin/xmgfas.out ]; then \
	../bin/xmgfas.out > xmgfas.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmgfas.out ../src/xmgfas.f $(ILIB) && \
	xmgfas.out > xmgfas.reslt && /bin/rm -f xmgfas.o xmgfas.out; \
	fi

xmglin.reslt : ../src/xmglin.f $(ILIB)
	-if [ -r ../bin/xmglin.out ]; then \
	../bin/xmglin.out > xmglin.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmglin.out ../src/xmglin.f $(ILIB) && \
	xmglin.out > xmglin.reslt && /bin/rm -f xmglin.o xmglin.out; \
	fi

xmidpnt.reslt : ../src/xmidpnt.f $(ILIB)
	-if [ -r ../bin/xmidpnt.out ]; then \
	../bin/xmidpnt.out > xmidpnt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmidpnt.out ../src/xmidpnt.f $(ILIB) && \
	xmidpnt.out > xmidpnt.reslt && /bin/rm -f xmidpnt.o xmidpnt.out; \
	fi

xmiser.reslt : ../src/xmiser.f $(ILIB)
	-if [ -r ../bin/xmiser.out ]; then \
	../bin/xmiser.out < ../responses/xmiser.resp > xmiser.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmiser.out ../src/xmiser.f $(ILIB) && \
	xmiser.out < ../responses/xmiser.resp > xmiser.reslt && \
	/bin/rm -f xmiser.o xmiser.out; \
	fi

xmmid.reslt : ../src/xmmid.f $(ILIB)
	-if [ -r ../bin/xmmid.out ]; then \
	../bin/xmmid.out < ../responses/xmmid.resp > xmmid.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmmid.out ../src/xmmid.f $(ILIB) && \
	xmmid.out < ../responses/xmmid.resp > xmmid.reslt && \
	/bin/rm -f xmmid.o xmmid.out; \
	fi

xmnbrak.reslt : ../src/xmnbrak.f $(ILIB)
	-if [ -r ../bin/xmnbrak.out ]; then \
	../bin/xmnbrak.out > xmnbrak.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmnbrak.out ../src/xmnbrak.f $(ILIB) && \
	xmnbrak.out > xmnbrak.reslt && /bin/rm -f xmnbrak.o xmnbrak.out; \
	fi

xmnewt.reslt : ../src/xmnewt.f $(ILIB)
	-if [ -r ../bin/xmnewt.out ]; then \
	../bin/xmnewt.out < ../responses/xmnewt.resp > xmnewt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmnewt.out ../src/xmnewt.f $(ILIB) && \
	xmnewt.out < ../responses/xmnewt.resp > xmnewt.reslt && \
	/bin/rm -f xmnewt.o xmnewt.out; \
	fi

xmoment.reslt : ../src/xmoment.f $(ILIB)
	-if [ -r ../bin/xmoment.out ]; then \
	../bin/xmoment.out > xmoment.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmoment.out ../src/xmoment.f $(ILIB) && \
	xmoment.out > xmoment.reslt && /bin/rm -f xmoment.o xmoment.out; \
	fi

xmppi.reslt : ../src/xmppi.f $(ILIB)
	-if [ -r ../bin/xmppi.out ]; then \
	../bin/xmppi.out < ../responses/xmppi.resp > xmppi.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmppi.out ../src/xmppi.f $(ILIB) && \
	xmppi.out < ../responses/xmppi.resp > xmppi.reslt && \
	/bin/rm -f xmppi.o xmppi.out; \
	fi

xmprove.reslt : ../src/xmprove.f $(ILIB)
	-if [ -r ../bin/xmprove.out ]; then \
	../bin/xmprove.out > xmprove.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmprove.out ../src/xmprove.f $(ILIB) && \
	xmprove.out > xmprove.reslt && /bin/rm -f xmprove.o xmprove.out; \
	fi

xmrqcof.reslt : ../src/xmrqcof.f $(ILIB)
	-if [ -r ../bin/xmrqcof.out ]; then \
	../bin/xmrqcof.out > xmrqcof.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmrqcof.out ../src/xmrqcof.f $(ILIB) && \
	xmrqcof.out > xmrqcof.reslt && /bin/rm -f xmrqcof.o xmrqcof.out; \
	fi

xmrqmin.reslt : ../src/xmrqmin.f $(ILIB)
	-if [ -r ../bin/xmrqmin.out ]; then \
	../bin/xmrqmin.out < ../responses/xmrqmin.resp > xmrqmin.reslt; \
	else \
	$(FC) $(FFLAGS) -o xmrqmin.out ../src/xmrqmin.f $(ILIB) && \
	xmrqmin.out < ../responses/xmrqmin.resp > xmrqmin.reslt && \
	/bin/rm -f xmrqmin.o xmrqmin.out; \
	fi

xnewt.reslt : ../src/xnewt.f $(ILIB)
	-if [ -r ../bin/xnewt.out ]; then \
	../bin/xnewt.out > xnewt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xnewt.out ../src/xnewt.f $(ILIB) && \
	xnewt.out > xnewt.reslt && /bin/rm -f xnewt.o xnewt.out; \
	fi

xodeint.reslt : ../src/xodeint.f $(ILIB)
	-if [ -r ../bin/xodeint.out ]; then \
	../bin/xodeint.out > xodeint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xodeint.out ../src/xodeint.f $(ILIB) && \
	xodeint.out > xodeint.reslt && /bin/rm -f xodeint.o xodeint.out; \
	fi

xorthog.reslt : ../src/xorthog.f $(ILIB)
	-if [ -r ../bin/xorthog.out ]; then \
	../bin/xorthog.out < ../responses/xorthog.resp > xorthog.reslt; \
	else \
	$(FC) $(FFLAGS) -o xorthog.out ../src/xorthog.f $(ILIB) && \
	xorthog.out < ../responses/xorthog.resp > xorthog.reslt && \
	/bin/rm -f xorthog.o xorthog.out; \
	fi

xpade.reslt : ../src/xpade.f $(ILIB)
	-if [ -r ../bin/xpade.out ]; then \
	../bin/xpade.out < ../responses/xpade.resp > xpade.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpade.out ../src/xpade.f $(ILIB) && \
	xpade.out < ../responses/xpade.resp > xpade.reslt && \
	/bin/rm -f xpade.o xpade.out; \
	fi

xpccheb.reslt : ../src/xpccheb.f $(ILIB)
	-if [ -r ../bin/xpccheb.out ]; then \
	../bin/xpccheb.out > xpccheb.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpccheb.out ../src/xpccheb.f $(ILIB) && \
	xpccheb.out > xpccheb.reslt && /bin/rm -f xpccheb.o xpccheb.out; \
	fi

xpcshft.reslt : ../src/xpcshft.f $(ILIB)
	-if [ -r ../bin/xpcshft.out ]; then \
	../bin/xpcshft.out < ../responses/xpcshft.resp > xpcshft.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpcshft.out ../src/xpcshft.f $(ILIB) && \
	xpcshft.out < ../responses/xpcshft.resp > xpcshft.reslt && \
	/bin/rm -f xpcshft.o xpcshft.out; \
	fi

xpearsn.reslt : ../src/xpearsn.f $(ILIB)
	-if [ -r ../bin/xpearsn.out ]; then \
	../bin/xpearsn.out > xpearsn.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpearsn.out ../src/xpearsn.f $(ILIB) && \
	xpearsn.out > xpearsn.reslt && /bin/rm -f xpearsn.o xpearsn.out; \
	fi

xperiod.reslt : ../src/xperiod.f $(ILIB)
	-if [ -r ../bin/xperiod.out ]; then \
	../bin/xperiod.out > xperiod.reslt; \
	else \
	$(FC) $(FFLAGS) -o xperiod.out ../src/xperiod.f $(ILIB) && \
	xperiod.out > xperiod.reslt && /bin/rm -f xperiod.o xperiod.out; \
	fi

xpiksr2.reslt : ../src/xpiksr2.f $(ILIB)
	-if [ -r ../bin/xpiksr2.out ]; then \
	../bin/xpiksr2.out < ../responses/xpiksr2.resp > xpiksr2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpiksr2.out ../src/xpiksr2.f $(ILIB) && \
	xpiksr2.out < ../responses/xpiksr2.resp > xpiksr2.reslt && \
	/bin/rm -f xpiksr2.o xpiksr2.out; \
	fi

xpiksrt.reslt : ../src/xpiksrt.f $(ILIB)
	-if [ -r ../bin/xpiksrt.out ]; then \
	../bin/xpiksrt.out > xpiksrt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpiksrt.out ../src/xpiksrt.f $(ILIB) && \
	xpiksrt.out > xpiksrt.reslt && /bin/rm -f xpiksrt.o xpiksrt.out; \
	fi

xplgndr.reslt : ../src/xplgndr.f $(ILIB)
	-if [ -r ../bin/xplgndr.out ]; then \
	../bin/xplgndr.out > xplgndr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xplgndr.out ../src/xplgndr.f $(ILIB) && \
	xplgndr.out > xplgndr.reslt && /bin/rm -f xplgndr.o xplgndr.out; \
	fi

xpoidev.reslt : ../src/xpoidev.f $(ILIB)
	-if [ -r ../bin/xpoidev.out ]; then \
	../bin/xpoidev.out < ../responses/xpoidev.resp > xpoidev.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpoidev.out ../src/xpoidev.f $(ILIB) && \
	xpoidev.out < ../responses/xpoidev.resp > xpoidev.reslt && \
	/bin/rm -f xpoidev.o xpoidev.out; \
	fi

xpolcoe.reslt : ../src/xpolcoe.f $(ILIB)
	-if [ -r ../bin/xpolcoe.out ]; then \
	../bin/xpolcoe.out < ../responses/xpolcoe.resp > xpolcoe.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpolcoe.out ../src/xpolcoe.f $(ILIB) && \
	xpolcoe.out < ../responses/xpolcoe.resp > xpolcoe.reslt && \
	/bin/rm -f xpolcoe.o xpolcoe.out; \
	fi

xpolcof.reslt : ../src/xpolcof.f $(ILIB)
	-if [ -r ../bin/xpolcof.out ]; then \
	../bin/xpolcof.out < ../responses/xpolcof.resp > xpolcof.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpolcof.out ../src/xpolcof.f $(ILIB) && \
	xpolcof.out < ../responses/xpolcof.resp > xpolcof.reslt && \
	/bin/rm -f xpolcof.o xpolcof.out; \
	fi

xpoldiv.reslt : ../src/xpoldiv.f $(ILIB)
	-if [ -r ../bin/xpoldiv.out ]; then \
	../bin/xpoldiv.out > xpoldiv.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpoldiv.out ../src/xpoldiv.f $(ILIB) && \
	xpoldiv.out > xpoldiv.reslt && /bin/rm -f xpoldiv.o xpoldiv.out; \
	fi

xpolin2.reslt : ../src/xpolin2.f $(ILIB)
	-if [ -r ../bin/xpolin2.out ]; then \
	../bin/xpolin2.out > xpolin2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpolin2.out ../src/xpolin2.f $(ILIB) && \
	xpolin2.out > xpolin2.reslt && /bin/rm -f xpolin2.o xpolin2.out; \
	fi

xpolint.reslt : ../src/xpolint.f $(ILIB)
	-if [ -r ../bin/xpolint.out ]; then \
	../bin/xpolint.out < ../responses/xpolint.resp > xpolint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpolint.out ../src/xpolint.f $(ILIB) && \
	xpolint.out < ../responses/xpolint.resp > xpolint.reslt && \
	/bin/rm -f xpolint.o xpolint.out; \
	fi

xpowell.reslt : ../src/xpowell.f $(ILIB)
	-if [ -r ../bin/xpowell.out ]; then \
	../bin/xpowell.out > xpowell.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpowell.out ../src/xpowell.f $(ILIB) && \
	xpowell.out > xpowell.reslt && /bin/rm -f xpowell.o xpowell.out; \
	fi

xpredic.reslt : ../src/xpredic.f $(ILIB)
	-if [ -r ../bin/xpredic.out ]; then \
	../bin/xpredic.out > xpredic.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpredic.out ../src/xpredic.f $(ILIB) && \
	xpredic.out > xpredic.reslt && /bin/rm -f xpredic.o xpredic.out; \
	fi

xprobks.reslt : ../src/xprobks.f $(ILIB)
	-if [ -r ../bin/xprobks.out ]; then \
	../bin/xprobks.out > xprobks.reslt; \
	else \
	$(FC) $(FFLAGS) -o xprobks.out ../src/xprobks.f $(ILIB) && \
	xprobks.out > xprobks.reslt && /bin/rm -f xprobks.o xprobks.out; \
	fi

xpsdes.reslt : ../src/xpsdes.f $(ILIB)
	-if [ -r ../bin/xpsdes.out ]; then \
	../bin/xpsdes.out > xpsdes.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpsdes.out ../src/xpsdes.f $(ILIB) && \
	xpsdes.out > xpsdes.reslt && /bin/rm -f xpsdes.o xpsdes.out; \
	fi

xpzextr.reslt : ../src/xpzextr.f $(ILIB)
	-if [ -r ../bin/xpzextr.out ]; then \
	../bin/xpzextr.out > xpzextr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xpzextr.out ../src/xpzextr.f $(ILIB) && \
	xpzextr.out > xpzextr.reslt && /bin/rm -f xpzextr.o xpzextr.out; \
	fi

xqgaus.reslt : ../src/xqgaus.f $(ILIB)
	-if [ -r ../bin/xqgaus.out ]; then \
	../bin/xqgaus.out > xqgaus.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqgaus.out ../src/xqgaus.f $(ILIB) && \
	xqgaus.out > xqgaus.reslt && /bin/rm -f xqgaus.o xqgaus.out; \
	fi

xqrdcmp.reslt : ../src/xqrdcmp.f $(ILIB)
	-if [ -r ../bin/xqrdcmp.out ]; then \
	../bin/xqrdcmp.out < ../responses/xqrdcmp.resp > xqrdcmp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqrdcmp.out ../src/xqrdcmp.f $(ILIB) && \
	xqrdcmp.out < ../responses/xqrdcmp.resp > xqrdcmp.reslt && \
	/bin/rm -f xqrdcmp.o xqrdcmp.out; \
	fi

xqromb.reslt : ../src/xqromb.f $(ILIB)
	-if [ -r ../bin/xqromb.out ]; then \
	../bin/xqromb.out > xqromb.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqromb.out ../src/xqromb.f $(ILIB) && \
	xqromb.out > xqromb.reslt && /bin/rm -f xqromb.o xqromb.out; \
	fi

xqromo.reslt : ../src/xqromo.f $(ILIB)
	-if [ -r ../bin/xqromo.out ]; then \
	../bin/xqromo.out > xqromo.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqromo.out ../src/xqromo.f $(ILIB) && \
	xqromo.out > xqromo.reslt && /bin/rm -f xqromo.o xqromo.out; \
	fi

xqroot.reslt : ../src/xqroot.f $(ILIB)
	-if [ -r ../bin/xqroot.out ]; then \
	../bin/xqroot.out > xqroot.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqroot.out ../src/xqroot.f $(ILIB) && \
	xqroot.out > xqroot.reslt && /bin/rm -f xqroot.o xqroot.out; \
	fi

xqrsolv.reslt : ../src/xqrsolv.f $(ILIB)
	-if [ -r ../bin/xqrsolv.out ]; then \
	../bin/xqrsolv.out < ../responses/xqrsolv.resp > xqrsolv.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqrsolv.out ../src/xqrsolv.f $(ILIB) && \
	xqrsolv.out < ../responses/xqrsolv.resp > xqrsolv.reslt && \
	/bin/rm -f xqrsolv.o xqrsolv.out; \
	fi

xqrupdt.reslt : ../src/xqrupdt.f $(ILIB)
	-if [ -r ../bin/xqrupdt.out ]; then \
	../bin/xqrupdt.out < ../responses/xqrupdt.resp > xqrupdt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqrupdt.out ../src/xqrupdt.f $(ILIB) && \
	xqrupdt.out < ../responses/xqrupdt.resp > xqrupdt.reslt && \
	/bin/rm -f xqrupdt.o xqrupdt.out; \
	fi

xqsimp.reslt : ../src/xqsimp.f $(ILIB)
	-if [ -r ../bin/xqsimp.out ]; then \
	../bin/xqsimp.out > xqsimp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqsimp.out ../src/xqsimp.f $(ILIB) && \
	xqsimp.out > xqsimp.reslt && /bin/rm -f xqsimp.o xqsimp.out; \
	fi

xqtrap.reslt : ../src/xqtrap.f $(ILIB)
	-if [ -r ../bin/xqtrap.out ]; then \
	../bin/xqtrap.out > xqtrap.reslt; \
	else \
	$(FC) $(FFLAGS) -o xqtrap.out ../src/xqtrap.f $(ILIB) && \
	xqtrap.out > xqtrap.reslt && /bin/rm -f xqtrap.o xqtrap.out; \
	fi

xquad3d.reslt : ../src/xquad3d.f $(ILIB)
	-if [ -r ../bin/xquad3d.out ]; then \
	../bin/xquad3d.out > xquad3d.reslt; \
	else \
	$(FC) $(FFLAGS) -o xquad3d.out ../src/xquad3d.f $(ILIB) && \
	xquad3d.out > xquad3d.reslt && /bin/rm -f xquad3d.o xquad3d.out; \
	fi

xran.reslt : ../src/xran.f $(ILIB)
	-if [ -r ../bin/xran.out ]; then \
	../bin/xran.out > xran.reslt; \
	else \
	$(FC) $(FFLAGS) -o xran.out ../src/xran.f $(ILIB) && \
	xran.out > xran.reslt && /bin/rm -f xran.o xran.out; \
	fi

xran4.reslt : ../src/xran4.f $(ILIB)
	-if [ -r ../bin/xran4.out ]; then \
	../bin/xran4.out > xran4.reslt; \
	else \
	$(FC) $(FFLAGS) -o xran4.out ../src/xran4.f $(ILIB) && \
	xran4.out > xran4.reslt && /bin/rm -f xran4.o xran4.out; \
	fi

xrank.reslt : ../src/xrank.f $(ILIB)
	-if [ -r ../bin/xrank.out ]; then \
	../bin/xrank.out < ../responses/xrank.resp > xrank.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrank.out ../src/xrank.f $(ILIB) && \
	xrank.out < ../responses/xrank.resp > xrank.reslt && \
	/bin/rm -f xrank.o xrank.out; \
	fi

xratint.reslt : ../src/xratint.f $(ILIB)
	-if [ -r ../bin/xratint.out ]; then \
	../bin/xratint.out > xratint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xratint.out ../src/xratint.f $(ILIB) && \
	xratint.out > xratint.reslt && /bin/rm -f xratint.o xratint.out; \
	fi

xratlsq.reslt : ../src/xratlsq.f $(ILIB)
	-if [ -r ../bin/xratlsq.out ]; then \
	../bin/xratlsq.out < ../responses/xratlsq.resp > xratlsq.reslt; \
	else \
	$(FC) $(FFLAGS) -o xratlsq.out ../src/xratlsq.f $(ILIB) && \
	xratlsq.out < ../responses/xratlsq.resp > xratlsq.reslt && \
	/bin/rm -f xratlsq.o xratlsq.out; \
	fi

xrc.reslt : ../src/xrc.f $(ILIB)
	-if [ -r ../bin/xrc.out ]; then \
	../bin/xrc.out > xrc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrc.out ../src/xrc.f $(ILIB) && \
	xrc.out > xrc.reslt && /bin/rm -f xrc.o xrc.out; \
	fi

xrd.reslt : ../src/xrd.f $(ILIB)
	-if [ -r ../bin/xrd.out ]; then \
	../bin/xrd.out > xrd.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrd.out ../src/xrd.f $(ILIB) && \
	xrd.out > xrd.reslt && /bin/rm -f xrd.o xrd.out; \
	fi

xrealft.reslt : ../src/xrealft.f $(ILIB)
	-if [ -r ../bin/xrealft.out ]; then \
	../bin/xrealft.out < ../responses/xrealft.resp > xrealft.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrealft.out ../src/xrealft.f $(ILIB) && \
	xrealft.out < ../responses/xrealft.resp > xrealft.reslt && \
	/bin/rm -f xrealft.o xrealft.out; \
	fi

xrf.reslt : ../src/xrf.f $(ILIB)
	-if [ -r ../bin/xrf.out ]; then \
	../bin/xrf.out > xrf.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrf.out ../src/xrf.f $(ILIB) && \
	xrf.out > xrf.reslt && /bin/rm -f xrf.o xrf.out; \
	fi

xrj.reslt : ../src/xrj.f $(ILIB)
	-if [ -r ../bin/xrj.out ]; then \
	../bin/xrj.out > xrj.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrj.out ../src/xrj.f $(ILIB) && \
	xrj.out > xrj.reslt && /bin/rm -f xrj.o xrj.out; \
	fi

xrk4.reslt : ../src/xrk4.f $(ILIB)
	-if [ -r ../bin/xrk4.out ]; then \
	../bin/xrk4.out > xrk4.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrk4.out ../src/xrk4.f $(ILIB) && \
	xrk4.out > xrk4.reslt && /bin/rm -f xrk4.o xrk4.out; \
	fi

xrkdumb.reslt : ../src/xrkdumb.f $(ILIB)
	-if [ -r ../bin/xrkdumb.out ]; then \
	../bin/xrkdumb.out > xrkdumb.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrkdumb.out ../src/xrkdumb.f $(ILIB) && \
	xrkdumb.out > xrkdumb.reslt && /bin/rm -f xrkdumb.o xrkdumb.out; \
	fi

xrkqs.reslt : ../src/xrkqs.f $(ILIB)
	-if [ -r ../bin/xrkqs.out ]; then \
	../bin/xrkqs.out > xrkqs.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrkqs.out ../src/xrkqs.f $(ILIB) && \
	xrkqs.out > xrkqs.reslt && /bin/rm -f xrkqs.o xrkqs.out; \
	fi

xrlft3.reslt : ../src/xrlft3.f $(ILIB)
	-if [ -r ../bin/xrlft3.out ]; then \
	../bin/xrlft3.out > xrlft3.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrlft3.out ../src/xrlft3.f $(ILIB) && \
	xrlft3.out > xrlft3.reslt && /bin/rm -f xrlft3.o xrlft3.out; \
	fi

xrofunc.reslt : ../src/xrofunc.f $(ILIB)
	-if [ -r ../bin/xrofunc.out ]; then \
	../bin/xrofunc.out > xrofunc.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrofunc.out ../src/xrofunc.f $(ILIB) && \
	xrofunc.out > xrofunc.reslt && /bin/rm -f xrofunc.o xrofunc.out; \
	fi

xrtbis.reslt : ../src/xrtbis.f $(ILIB)
	-if [ -r ../bin/xrtbis.out ]; then \
	../bin/xrtbis.out > xrtbis.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrtbis.out ../src/xrtbis.f $(ILIB) && \
	xrtbis.out > xrtbis.reslt && /bin/rm -f xrtbis.o xrtbis.out; \
	fi

xrtflsp.reslt : ../src/xrtflsp.f $(ILIB)
	-if [ -r ../bin/xrtflsp.out ]; then \
	../bin/xrtflsp.out > xrtflsp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrtflsp.out ../src/xrtflsp.f $(ILIB) && \
	xrtflsp.out > xrtflsp.reslt && /bin/rm -f xrtflsp.o xrtflsp.out; \
	fi

xrtnewt.reslt : ../src/xrtnewt.f $(ILIB)
	-if [ -r ../bin/xrtnewt.out ]; then \
	../bin/xrtnewt.out > xrtnewt.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrtnewt.out ../src/xrtnewt.f $(ILIB) && \
	xrtnewt.out > xrtnewt.reslt && /bin/rm -f xrtnewt.o xrtnewt.out; \
	fi

xrtsafe.reslt : ../src/xrtsafe.f $(ILIB)
	-if [ -r ../bin/xrtsafe.out ]; then \
	../bin/xrtsafe.out > xrtsafe.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrtsafe.out ../src/xrtsafe.f $(ILIB) && \
	xrtsafe.out > xrtsafe.reslt && /bin/rm -f xrtsafe.o xrtsafe.out; \
	fi

xrtsec.reslt : ../src/xrtsec.f $(ILIB)
	-if [ -r ../bin/xrtsec.out ]; then \
	../bin/xrtsec.out > xrtsec.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrtsec.out ../src/xrtsec.f $(ILIB) && \
	xrtsec.out > xrtsec.reslt && /bin/rm -f xrtsec.o xrtsec.out; \
	fi

xrzextr.reslt : ../src/xrzextr.f $(ILIB)
	-if [ -r ../bin/xrzextr.out ]; then \
	../bin/xrzextr.out > xrzextr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xrzextr.out ../src/xrzextr.f $(ILIB) && \
	xrzextr.out > xrzextr.reslt && /bin/rm -f xrzextr.o xrzextr.out; \
	fi

xsavgol.reslt : ../src/xsavgol.f $(ILIB)
	-if [ -r ../bin/xsavgol.out ]; then \
	../bin/xsavgol.out > xsavgol.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsavgol.out ../src/xsavgol.f $(ILIB) && \
	xsavgol.out > xsavgol.reslt && /bin/rm -f xsavgol.o xsavgol.out; \
	fi

xscrsho.reslt : ../src/xscrsho.f $(ILIB)
	-if [ -r ../bin/xscrsho.out ]; then \
	../bin/xscrsho.out < ../responses/xscrsho.resp > xscrsho.reslt; \
	else \
	$(FC) $(FFLAGS) -o xscrsho.out ../src/xscrsho.f $(ILIB) && \
	xscrsho.out < ../responses/xscrsho.resp > xscrsho.reslt && \
	/bin/rm -f xscrsho.o xscrsho.out; \
	fi

xselect.reslt : ../src/xselect.f $(ILIB)
	-if [ -r ../bin/xselect.out ]; then \
	../bin/xselect.out < ../responses/xselect.resp > xselect.reslt; \
	else \
	$(FC) $(FFLAGS) -o xselect.out ../src/xselect.f $(ILIB) && \
	xselect.out < ../responses/xselect.resp > xselect.reslt && \
	/bin/rm -f xselect.o xselect.out; \
	fi

xselip.reslt : ../src/xselip.f $(ILIB)
	-if [ -r ../bin/xselip.out ]; then \
	../bin/xselip.out < ../responses/xselip.resp > xselip.reslt; \
	else \
	$(FC) $(FFLAGS) -o xselip.out ../src/xselip.f $(ILIB) && \
	xselip.out < ../responses/xselip.resp > xselip.reslt && \
	/bin/rm -f xselip.o xselip.out; \
	fi

xshell.reslt : ../src/xshell.f $(ILIB)
	-if [ -r ../bin/xshell.out ]; then \
	../bin/xshell.out > xshell.reslt; \
	else \
	$(FC) $(FFLAGS) -o xshell.out ../src/xshell.f $(ILIB) && \
	xshell.out > xshell.reslt && /bin/rm -f xshell.o xshell.out; \
	fi

xsimplx.reslt : ../src/xsimplx.f $(ILIB)
	-if [ -r ../bin/xsimplx.out ]; then \
	../bin/xsimplx.out > xsimplx.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsimplx.out ../src/xsimplx.f $(ILIB) && \
	xsimplx.out > xsimplx.reslt && /bin/rm -f xsimplx.o xsimplx.out; \
	fi

xsimpr.reslt : ../src/xsimpr.f $(ILIB)
	-if [ -r ../bin/xsimpr.out ]; then \
	../bin/xsimpr.out > xsimpr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsimpr.out ../src/xsimpr.f $(ILIB) && \
	xsimpr.out > xsimpr.reslt && /bin/rm -f xsimpr.o xsimpr.out; \
	fi

xsinft.reslt : ../src/xsinft.f $(ILIB)
	-if [ -r ../bin/xsinft.out ]; then \
	../bin/xsinft.out < ../responses/xsinft.resp > xsinft.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsinft.out ../src/xsinft.f $(ILIB) && \
	xsinft.out < ../responses/xsinft.resp > xsinft.reslt && \
	/bin/rm -f xsinft.o xsinft.out; \
	fi

xsncndn.reslt : ../src/xsncndn.f $(ILIB)
	-if [ -r ../bin/xsncndn.out ]; then \
	../bin/xsncndn.out > xsncndn.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsncndn.out ../src/xsncndn.f $(ILIB) && \
	xsncndn.out > xsncndn.reslt && /bin/rm -f xsncndn.o xsncndn.out; \
	fi

xsobseq.reslt : ../src/xsobseq.f $(ILIB)
	-if [ -r ../bin/xsobseq.out ]; then \
	../bin/xsobseq.out > xsobseq.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsobseq.out ../src/xsobseq.f $(ILIB) && \
	xsobseq.out > xsobseq.reslt && /bin/rm -f xsobseq.o xsobseq.out; \
	fi

xsor.reslt : ../src/xsor.f $(ILIB)
	-if [ -r ../bin/xsor.out ]; then \
	../bin/xsor.out > xsor.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsor.out ../src/xsor.f $(ILIB) && \
	xsor.out > xsor.reslt && /bin/rm -f xsor.o xsor.out; \
	fi

xsort.reslt : ../src/xsort.f $(ILIB)
	-if [ -r ../bin/xsort.out ]; then \
	../bin/xsort.out > xsort.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsort.out ../src/xsort.f $(ILIB) && \
	xsort.out > xsort.reslt && /bin/rm -f xsort.o xsort.out; \
	fi

xsort2.reslt : ../src/xsort2.f $(ILIB)
	-if [ -r ../bin/xsort2.out ]; then \
	../bin/xsort2.out < ../responses/xsort2.resp > xsort2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsort2.out ../src/xsort2.f $(ILIB) && \
	xsort2.out < ../responses/xsort2.resp > xsort2.reslt && \
	/bin/rm -f xsort2.o xsort2.out; \
	fi

xsort3.reslt : ../src/xsort3.f $(ILIB)
	-if [ -r ../bin/xsort3.out ]; then \
	../bin/xsort3.out > xsort3.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsort3.out ../src/xsort3.f $(ILIB) && \
	xsort3.out > xsort3.reslt && /bin/rm -f xsort3.o xsort3.out; \
	fi

xspctrm.reslt : ../src/xspctrm.f $(ILIB)
	-if [ -r ../bin/xspctrm.out ]; then \
	../bin/xspctrm.out > xspctrm.reslt; \
	else \
	$(FC) $(FFLAGS) -o xspctrm.out ../src/xspctrm.f $(ILIB) && \
	xspctrm.out > xspctrm.reslt && /bin/rm -f xspctrm.o xspctrm.out; \
	fi

xspear.reslt : ../src/xspear.f $(ILIB)
	-if [ -r ../bin/xspear.out ]; then \
	../bin/xspear.out > xspear.reslt; \
	else \
	$(FC) $(FFLAGS) -o xspear.out ../src/xspear.f $(ILIB) && \
	xspear.out > xspear.reslt && /bin/rm -f xspear.o xspear.out; \
	fi

xsphbes.reslt : ../src/xsphbes.f $(ILIB)
	-if [ -r ../bin/xsphbes.out ]; then \
	../bin/xsphbes.out > xsphbes.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsphbes.out ../src/xsphbes.f $(ILIB) && \
	xsphbes.out > xsphbes.reslt && /bin/rm -f xsphbes.o xsphbes.out; \
	fi

xsplie2.reslt : ../src/xsplie2.f $(ILIB)
	-if [ -r ../bin/xsplie2.out ]; then \
	../bin/xsplie2.out > xsplie2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsplie2.out ../src/xsplie2.f $(ILIB) && \
	xsplie2.out > xsplie2.reslt && /bin/rm -f xsplie2.o xsplie2.out; \
	fi

xsplin2.reslt : ../src/xsplin2.f $(ILIB)
	-if [ -r ../bin/xsplin2.out ]; then \
	../bin/xsplin2.out > xsplin2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsplin2.out ../src/xsplin2.f $(ILIB) && \
	xsplin2.out > xsplin2.reslt && /bin/rm -f xsplin2.o xsplin2.out; \
	fi

xspline.reslt : ../src/xspline.f $(ILIB)
	-if [ -r ../bin/xspline.out ]; then \
	../bin/xspline.out > xspline.reslt; \
	else \
	$(FC) $(FFLAGS) -o xspline.out ../src/xspline.f $(ILIB) && \
	xspline.out > xspline.reslt && /bin/rm -f xspline.o xspline.out; \
	fi

xsplint.reslt : ../src/xsplint.f $(ILIB)
	-if [ -r ../bin/xsplint.out ]; then \
	../bin/xsplint.out < ../responses/xsplint.resp > xsplint.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsplint.out ../src/xsplint.f $(ILIB) && \
	xsplint.out < ../responses/xsplint.resp > xsplint.reslt && \
	/bin/rm -f xsplint.o xsplint.out; \
	fi

xsprsax.reslt : ../src/xsprsax.f $(ILIB)
	-if [ -r ../bin/xsprsax.out ]; then \
	../bin/xsprsax.out > xsprsax.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsprsax.out ../src/xsprsax.f $(ILIB) && \
	xsprsax.out > xsprsax.reslt && /bin/rm -f xsprsax.o xsprsax.out; \
	fi

xsprsin.reslt : ../src/xsprsin.f $(ILIB)
	-if [ -r ../bin/xsprsin.out ]; then \
	../bin/xsprsin.out > xsprsin.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsprsin.out ../src/xsprsin.f $(ILIB) && \
	xsprsin.out > xsprsin.reslt && /bin/rm -f xsprsin.o xsprsin.out; \
	fi

xsprspm.reslt : ../src/xsprspm.f $(ILIB)
	-if [ -r ../bin/xsprspm.out ]; then \
	../bin/xsprspm.out > xsprspm.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsprspm.out ../src/xsprspm.f $(ILIB) && \
	xsprspm.out > xsprspm.reslt && /bin/rm -f xsprspm.o xsprspm.out; \
	fi

xsprstm.reslt : ../src/xsprstm.f $(ILIB)
	-if [ -r ../bin/xsprstm.out ]; then \
	../bin/xsprstm.out > xsprstm.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsprstm.out ../src/xsprstm.f $(ILIB) && \
	xsprstm.out > xsprstm.reslt && /bin/rm -f xsprstm.o xsprstm.out; \
	fi

xsprstp.reslt : ../src/xsprstp.f $(ILIB)
	-if [ -r ../bin/xsprstp.out ]; then \
	../bin/xsprstp.out > xsprstp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsprstp.out ../src/xsprstp.f $(ILIB) && \
	xsprstp.out > xsprstp.reslt && /bin/rm -f xsprstp.o xsprstp.out; \
	fi

xsprstx.reslt : ../src/xsprstx.f $(ILIB)
	-if [ -r ../bin/xsprstx.out ]; then \
	../bin/xsprstx.out > xsprstx.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsprstx.out ../src/xsprstx.f $(ILIB) && \
	xsprstx.out > xsprstx.reslt && /bin/rm -f xsprstx.o xsprstx.out; \
	fi

xstifbs.reslt : ../src/xstifbs.f $(ILIB)
	-if [ -r ../bin/xstifbs.out ]; then \
	../bin/xstifbs.out < ../responses/xstifbs.resp > xstifbs.reslt; \
	else \
	$(FC) $(FFLAGS) -o xstifbs.out ../src/xstifbs.f $(ILIB) && \
	xstifbs.out < ../responses/xstifbs.resp > xstifbs.reslt && \
	/bin/rm -f xstifbs.o xstifbs.out; \
	fi

xstiff.reslt : ../src/xstiff.f $(ILIB)
	-if [ -r ../bin/xstiff.out ]; then \
	../bin/xstiff.out < ../responses/xstiff.resp > xstiff.reslt; \
	else \
	$(FC) $(FFLAGS) -o xstiff.out ../src/xstiff.f $(ILIB) && \
	xstiff.out < ../responses/xstiff.resp > xstiff.reslt && \
	/bin/rm -f xstiff.o xstiff.out; \
	fi

xstoerm.reslt : ../src/xstoerm.f $(ILIB)
	-if [ -r ../bin/xstoerm.out ]; then \
	../bin/xstoerm.out > xstoerm.reslt; \
	else \
	$(FC) $(FFLAGS) -o xstoerm.out ../src/xstoerm.f $(ILIB) && \
	xstoerm.out > xstoerm.reslt && /bin/rm -f xstoerm.o xstoerm.out; \
	fi

xsvbksb.reslt : ../src/xsvbksb.f $(ILIB)
	-if [ -r ../bin/xsvbksb.out ]; then \
	../bin/xsvbksb.out < ../responses/xsvbksb.resp > xsvbksb.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsvbksb.out ../src/xsvbksb.f $(ILIB) && \
	xsvbksb.out < ../responses/xsvbksb.resp > xsvbksb.reslt && \
	/bin/rm -f xsvbksb.o xsvbksb.out; \
	fi

xsvdcmp.reslt : ../src/xsvdcmp.f $(ILIB)
	-if [ -r ../bin/xsvdcmp.out ]; then \
	../bin/xsvdcmp.out < ../responses/xsvdcmp.resp > xsvdcmp.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsvdcmp.out ../src/xsvdcmp.f $(ILIB) && \
	xsvdcmp.out < ../responses/xsvdcmp.resp > xsvdcmp.reslt && \
	/bin/rm -f xsvdcmp.o xsvdcmp.out; \
	fi

xsvdfit.reslt : ../src/xsvdfit.f $(ILIB)
	-if [ -r ../bin/xsvdfit.out ]; then \
	../bin/xsvdfit.out > xsvdfit.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsvdfit.out ../src/xsvdfit.f $(ILIB) && \
	xsvdfit.out > xsvdfit.reslt && /bin/rm -f xsvdfit.o xsvdfit.out; \
	fi

xsvdvar.reslt : ../src/xsvdvar.f $(ILIB)
	-if [ -r ../bin/xsvdvar.out ]; then \
	../bin/xsvdvar.out > xsvdvar.reslt; \
	else \
	$(FC) $(FFLAGS) -o xsvdvar.out ../src/xsvdvar.f $(ILIB) && \
	xsvdvar.out > xsvdvar.reslt && /bin/rm -f xsvdvar.o xsvdvar.out; \
	fi

xtoeplz.reslt : ../src/xtoeplz.f $(ILIB)
	-if [ -r ../bin/xtoeplz.out ]; then \
	../bin/xtoeplz.out > xtoeplz.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtoeplz.out ../src/xtoeplz.f $(ILIB) && \
	xtoeplz.out > xtoeplz.reslt && /bin/rm -f xtoeplz.o xtoeplz.out; \
	fi

xtptest.reslt : ../src/xtptest.f $(ILIB)
	-if [ -r ../bin/xtptest.out ]; then \
	../bin/xtptest.out > xtptest.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtptest.out ../src/xtptest.f $(ILIB) && \
	xtptest.out > xtptest.reslt && /bin/rm -f xtptest.o xtptest.out; \
	fi

xtqli.reslt : ../src/xtqli.f $(ILIB)
	-if [ -r ../bin/xtqli.out ]; then \
	../bin/xtqli.out < ../responses/xtqli.resp > xtqli.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtqli.out ../src/xtqli.f $(ILIB) && \
	xtqli.out < ../responses/xtqli.resp > xtqli.reslt && \
	/bin/rm -f xtqli.o xtqli.out; \
	fi

xtrapzd.reslt : ../src/xtrapzd.f $(ILIB)
	-if [ -r ../bin/xtrapzd.out ]; then \
	../bin/xtrapzd.out > xtrapzd.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtrapzd.out ../src/xtrapzd.f $(ILIB) && \
	xtrapzd.out > xtrapzd.reslt && /bin/rm -f xtrapzd.o xtrapzd.out; \
	fi

xtred2.reslt : ../src/xtred2.f $(ILIB)
	-if [ -r ../bin/xtred2.out ]; then \
	../bin/xtred2.out > xtred2.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtred2.out ../src/xtred2.f $(ILIB) && \
	xtred2.out > xtred2.reslt && /bin/rm -f xtred2.o xtred2.out; \
	fi

xtridag.reslt : ../src/xtridag.f $(ILIB)
	-if [ -r ../bin/xtridag.out ]; then \
	../bin/xtridag.out < ../responses/xtridag.resp > xtridag.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtridag.out ../src/xtridag.f $(ILIB) && \
	xtridag.out < ../responses/xtridag.resp > xtridag.reslt && \
	/bin/rm -f xtridag.o xtridag.out; \
	fi

xttest.reslt : ../src/xttest.f $(ILIB)
	-if [ -r ../bin/xttest.out ]; then \
	../bin/xttest.out > xttest.reslt; \
	else \
	$(FC) $(FFLAGS) -o xttest.out ../src/xttest.f $(ILIB) && \
	xttest.out > xttest.reslt && /bin/rm -f xttest.o xttest.out; \
	fi

xtutest.reslt : ../src/xtutest.f $(ILIB)
	-if [ -r ../bin/xtutest.out ]; then \
	../bin/xtutest.out > xtutest.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtutest.out ../src/xtutest.f $(ILIB) && \
	xtutest.out > xtutest.reslt && /bin/rm -f xtutest.o xtutest.out; \
	fi

xtwofft.reslt : ../src/xtwofft.f $(ILIB)
	-if [ -r ../bin/xtwofft.out ]; then \
	../bin/xtwofft.out < ../responses/xtwofft.resp > xtwofft.reslt; \
	else \
	$(FC) $(FFLAGS) -o xtwofft.out ../src/xtwofft.f $(ILIB) && \
	xtwofft.out < ../responses/xtwofft.resp > xtwofft.reslt && \
	/bin/rm -f xtwofft.o xtwofft.out; \
	fi

xvander.reslt : ../src/xvander.f $(ILIB)
	-if [ -r ../bin/xvander.out ]; then \
	../bin/xvander.out > xvander.reslt; \
	else \
	$(FC) $(FFLAGS) -o xvander.out ../src/xvander.f $(ILIB) && \
	xvander.out > xvander.reslt && /bin/rm -f xvander.o xvander.out; \
	fi

xvegas.reslt : ../src/xvegas.f $(ILIB)
	-if [ -r ../bin/xvegas.out ]; then \
	../bin/xvegas.out < ../responses/xvegas.resp > xvegas.reslt; \
	else \
	$(FC) $(FFLAGS) -o xvegas.out ../src/xvegas.f $(ILIB) && \
	xvegas.out < ../responses/xvegas.resp > xvegas.reslt && \
	/bin/rm -f xvegas.o xvegas.out; \
	fi

xvoltra.reslt : ../src/xvoltra.f $(ILIB)
	-if [ -r ../bin/xvoltra.out ]; then \
	../bin/xvoltra.out > xvoltra.reslt; \
	else \
	$(FC) $(FFLAGS) -o xvoltra.out ../src/xvoltra.f $(ILIB) && \
	xvoltra.out > xvoltra.reslt && /bin/rm -f xvoltra.o xvoltra.out; \
	fi

xwt1.reslt : ../src/xwt1.f $(ILIB)
	-if [ -r ../bin/xwt1.out ]; then \
	../bin/xwt1.out < ../responses/xwt1.resp > xwt1.reslt; \
	else \
	$(FC) $(FFLAGS) -o xwt1.out ../src/xwt1.f $(ILIB) && \
	xwt1.out < ../responses/xwt1.resp > xwt1.reslt && \
	/bin/rm -f xwt1.o xwt1.out; \
	fi

xwtn.reslt : ../src/xwtn.f $(ILIB)
	-if [ -r ../bin/xwtn.out ]; then \
	../bin/xwtn.out > xwtn.reslt; \
	else \
	$(FC) $(FFLAGS) -o xwtn.out ../src/xwtn.f $(ILIB) && \
	xwtn.out > xwtn.reslt && /bin/rm -f xwtn.o xwtn.out; \
	fi

xzbrac.reslt : ../src/xzbrac.f $(ILIB)
	-if [ -r ../bin/xzbrac.out ]; then \
	../bin/xzbrac.out > xzbrac.reslt; \
	else \
	$(FC) $(FFLAGS) -o xzbrac.out ../src/xzbrac.f $(ILIB) && \
	xzbrac.out > xzbrac.reslt && /bin/rm -f xzbrac.o xzbrac.out; \
	fi

xzbrak.reslt : ../src/xzbrak.f $(ILIB)
	-if [ -r ../bin/xzbrak.out ]; then \
	../bin/xzbrak.out > xzbrak.reslt; \
	else \
	$(FC) $(FFLAGS) -o xzbrak.out ../src/xzbrak.f $(ILIB) && \
	xzbrak.out > xzbrak.reslt && /bin/rm -f xzbrak.o xzbrak.out; \
	fi

xzbrent.reslt : ../src/xzbrent.f $(ILIB)
	-if [ -r ../bin/xzbrent.out ]; then \
	../bin/xzbrent.out > xzbrent.reslt; \
	else \
	$(FC) $(FFLAGS) -o xzbrent.out ../src/xzbrent.f $(ILIB) && \
	xzbrent.out > xzbrent.reslt && /bin/rm -f xzbrent.o xzbrent.out; \
	fi

xzrhqr.reslt : ../src/xzrhqr.f $(ILIB)
	-if [ -r ../bin/xzrhqr.out ]; then \
	../bin/xzrhqr.out > xzrhqr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xzrhqr.out ../src/xzrhqr.f $(ILIB) && \
	xzrhqr.out > xzrhqr.reslt && /bin/rm -f xzrhqr.o xzrhqr.out; \
	fi

xzriddr.reslt : ../src/xzriddr.f $(ILIB)
	-if [ -r ../bin/xzriddr.out ]; then \
	../bin/xzriddr.out > xzriddr.reslt; \
	else \
	$(FC) $(FFLAGS) -o xzriddr.out ../src/xzriddr.f $(ILIB) && \
	xzriddr.out > xzriddr.reslt && /bin/rm -f xzriddr.o xzriddr.out; \
	fi

xzroots.reslt : ../src/xzroots.f $(ILIB)
	-if [ -r ../bin/xzroots.out ]; then \
	../bin/xzroots.out > xzroots.reslt; \
	else \
	$(FC) $(FFLAGS) -o xzroots.out ../src/xzroots.f $(ILIB) && \
	xzroots.out > xzroots.reslt && /bin/rm -f xzroots.o xzroots.out; \
	fi

badluk.reslt : ../../recipes/badluk.f $(ILIB)
	-if [ -r $(BINDIR)/badluk ]; then \
	$(BINDIR)/badluk > badluk.reslt; \
	else \
	$(FC) $(FFLAGS) -o badluk ../../recipes/badluk.f $(ILIB) && \
	badluk > badluk.reslt && /bin/rm -f badluk.o badluk; \
	fi

fredex.reslt : ../../recipes/fredex.f $(ILIB)
	-if [ -r $(BINDIR)/fredex ]; then \
	$(BINDIR)/fredex > fredex.reslt; \
	else \
	$(FC) $(FFLAGS) -o fredex ../../recipes/fredex.f $(ILIB) && \
	fredex > fredex.reslt && /bin/rm -f fredex.o fredex; \
	fi

sfroid.reslt : ../../recipes/sfroid.f $(ILIB)
	-if [ -r $(BINDIR)/sfroid ]; then \
	$(BINDIR)/sfroid < ../responses/sfroid.resp > sfroid.reslt; \
	else \
	$(FC) $(FFLAGS) -o sfroid ../../recipes/sfroid.f $(ILIB) && \
	sfroid < ../responses/sfroid.resp > sfroid.reslt && \
	/bin/rm -f sfroid.o sfroid; \
	fi

sphoot.reslt : ../../recipes/sphoot.f $(ILIB)
	-if [ -r $(BINDIR)/sphoot ]; then \
	$(BINDIR)/sphoot < ../responses/sphoot.resp > sphoot.reslt; \
	else \
	$(FC) $(FFLAGS) -o sphoot ../../recipes/sphoot.f \
	../../recipes/shoot.f $(ILIB) && \
	sphoot < ../responses/sphoot.resp > sphoot.reslt && \
	/bin/rm -f shoot.o sphoot.o sphoot; \
	fi

xsphfpt.reslt : ../../recipes/sphfpt.f ../src/xsphfpt.f $(ILIB)
	-if [ -r $(BINDIR)/sphfpt ]; then \
	$(BINDIR)/sphfpt < ../responses/xsphfpt.resp > xsphfpt.reslt; \
	else \
	$(FC) $(FFLAGS) -o sphfpt ../src/xsphfpt.f \
	../../recipes/sphfpt.f ../../recipes/shootf.f $(ILIB) && \
	sphfpt < ../responses/xsphfpt.resp > xsphfpt.reslt && \
	/bin/rm -f xsphfpt.o sphfpt.o shootf.o sphfpt; \
	fi
