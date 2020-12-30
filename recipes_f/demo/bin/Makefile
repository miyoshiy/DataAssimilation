# Makefile for Numerical Recipes Demo binaries

ILIB = $(LIBDIR)/$(FLIB)

OUTRULE = $(FC) $(FFLAGS) -o $@ $? $(ILIB)

all :	xairy.out xamebsa.out xamoeba.out xanneal.out xarcode.out \
	xavevar.out xbalanc.out xbandec.out xbanmul.out xbcucof.out \
	xbcuint.out xbeschb.out xbessi.out xbessi0.out xbessi1.out \
	xbessik.out xbessj.out xbessj0.out xbessj1.out xbessjy.out \
	xbessk.out xbessk0.out xbessk1.out xbessy.out xbessy0.out \
	xbessy1.out xbeta.out xbetai.out xbico.out xbnldev.out \
	xbrent.out xbroydn.out xbsstep.out xcaldat.out xchder.out \
	xchebev.out xchebft.out xchebpc.out xchint.out xcholsl.out \
	xchsone.out xchstwo.out xcisi.out xcntab1.out xcntab2.out \
	xconvlv.out xcorrel.out xcosft1.out xcosft2.out xcovsrt.out \
	xcrank.out xcyclic.out xdawson.out xdbrent.out xddpoly.out \
	xdecchk.out xdf1dim.out xdfpmin.out xdfridr.out xdftint.out \
	xeclass.out xeclazz.out xei.out xeigsrt.out xelle.out xellf.out \
	xellpi.out xelmhes.out xerf.out xerfc.out xerfcc.out xeulsum.out \
	xevlmem.out xexpdev.out xexpint.out xf1dim.out xfactln.out \
	xfactrl.out xfasper.out xfgauss.out xfit.out xfitexy.out \
	xfixrts.out xfleg.out xflmoon.out xfour1.out xfourfs.out \
	xfourn.out xfpoly.out xfred2.out xfredin.out xfrenel.out \
	xfrprmn.out xftest.out xgamdev.out xgammln.out xgammp.out \
	xgammq.out xgasdev.out xgaucof.out xgauher.out xgaujac.out \
	xgaulag.out xgauleg.out xgaussj.out xgcf.out xgolden.out \
	xgser.out xhpsel.out xhpsort.out xhqr.out xhuffman.out xhunt.out \
	xhypgeo.out xicrc.out xigray.out xindexx.out xirbit1.out \
	xirbit2.out xjacobi.out xjulday.out xkendl1.out xkendl2.out \
	xks2d1s.out xks2d2s.out xksone.out xkstwo.out xlaguer.out \
	xlfit.out xlinbcg.out xlinmin.out xlocate.out xlubksb.out \
	xludcmp.out xmachar.out xmedfit.out xmemcof.out xmgfas.out \
	xmglin.out xmidpnt.out xmiser.out xmmid.out xmnbrak.out \
	xmnewt.out xmoment.out xmppi.out xmprove.out xmrqcof.out \
	xmrqmin.out xnewt.out xodeint.out xorthog.out xpade.out \
	xpccheb.out xpcshft.out xpearsn.out xperiod.out xpiksr2.out \
	xpiksrt.out xplgndr.out xpoidev.out xpolcoe.out xpolcof.out \
	xpoldiv.out xpolin2.out xpolint.out xpowell.out xpredic.out \
	xprobks.out xpsdes.out xpzextr.out xqgaus.out xqrdcmp.out \
	xqromb.out xqromo.out xqroot.out xqrsolv.out xqrupdt.out \
	xqsimp.out xqtrap.out xquad3d.out xran.out xran4.out xrank.out \
	xratint.out xratlsq.out xrc.out xrd.out xrealft.out xrf.out \
	xrj.out xrk4.out xrkdumb.out xrkqs.out xrlft3.out xrofunc.out \
	xrtbis.out xrtflsp.out xrtnewt.out xrtsafe.out xrtsec.out \
	xrzextr.out xsavgol.out xscrsho.out xselect.out xselip.out \
	xshell.out xsimplx.out xsimpr.out xsinft.out xsncndn.out \
	xsobseq.out xsor.out xsort.out xsort2.out xsort3.out xspctrm.out \
	xspear.out xsphbes.out xsplie2.out xsplin2.out \
	xspline.out xsplint.out xsprsax.out xsprsin.out xsprspm.out \
	xsprstm.out xsprstp.out xsprstx.out xstifbs.out xstiff.out \
	xstoerm.out xsvbksb.out xsvdcmp.out xsvdfit.out xsvdvar.out \
	xtoeplz.out xtptest.out xtqli.out xtrapzd.out xtred2.out \
	xtridag.out xttest.out xtutest.out xtwofft.out xvander.out \
	xvegas.out xvoltra.out xwt1.out xwtn.out xzbrac.out xzbrak.out \
	xzbrent.out xzrhqr.out xzriddr.out xzroots.out

# svdcmp.dem.f doesn't compile right with optimization
xsvdcmp.out : ../src/xsvdcmp.f $(ILIB)
	case "$(FFLAGS)" in *-g*) $(FC) -g -o $@ $? $(ILIB) ;; \
	*) $(FC)  -o $@ $? $(ILIB) ;; esac
# If FFLAGS has -g in it, use the debug flag, else compile plain

xairy.out : ../src/xairy.f $(ILIB)
	$(OUTRULE)

xamebsa.out : ../src/xamebsa.f $(ILIB)
	$(OUTRULE)

xamoeba.out : ../src/xamoeba.f $(ILIB)
	$(OUTRULE)

xanneal.out : ../src/xanneal.f $(ILIB)
	$(OUTRULE)

xarcode.out : ../src/xarcode.f $(ILIB)
	$(OUTRULE)

xavevar.out : ../src/xavevar.f $(ILIB)
	$(OUTRULE)

xbalanc.out : ../src/xbalanc.f $(ILIB)
	$(OUTRULE)

xbandec.out : ../src/xbandec.f $(ILIB)
	$(OUTRULE)

xbanmul.out : ../src/xbanmul.f $(ILIB)
	$(OUTRULE)

xbcucof.out : ../src/xbcucof.f $(ILIB)
	$(OUTRULE)

xbcuint.out : ../src/xbcuint.f $(ILIB)
	$(OUTRULE)

xbeschb.out : ../src/xbeschb.f $(ILIB)
	$(OUTRULE)

xbessi.out : ../src/xbessi.f $(ILIB)
	$(OUTRULE)

xbessi0.out : ../src/xbessi0.f $(ILIB)
	$(OUTRULE)

xbessi1.out : ../src/xbessi1.f $(ILIB)
	$(OUTRULE)

xbessik.out : ../src/xbessik.f $(ILIB)
	$(OUTRULE)

xbessj.out : ../src/xbessj.f $(ILIB)
	$(OUTRULE)

xbessj0.out : ../src/xbessj0.f $(ILIB)
	$(OUTRULE)

xbessj1.out : ../src/xbessj1.f $(ILIB)
	$(OUTRULE)

xbessjy.out : ../src/xbessjy.f $(ILIB)
	$(OUTRULE)

xbessk.out : ../src/xbessk.f $(ILIB)
	$(OUTRULE)

xbessk0.out : ../src/xbessk0.f $(ILIB)
	$(OUTRULE)

xbessk1.out : ../src/xbessk1.f $(ILIB)
	$(OUTRULE)

xbessy.out : ../src/xbessy.f $(ILIB)
	$(OUTRULE)

xbessy0.out : ../src/xbessy0.f $(ILIB)
	$(OUTRULE)

xbessy1.out : ../src/xbessy1.f $(ILIB)
	$(OUTRULE)

xbeta.out : ../src/xbeta.f $(ILIB)
	$(OUTRULE)

xbetai.out : ../src/xbetai.f $(ILIB)
	$(OUTRULE)

xbico.out : ../src/xbico.f $(ILIB)
	$(OUTRULE)

xbnldev.out : ../src/xbnldev.f $(ILIB)
	$(OUTRULE)

xbrent.out : ../src/xbrent.f $(ILIB)
	$(OUTRULE)

xbroydn.out : ../src/xbroydn.f $(ILIB)
	$(OUTRULE)

xbsstep.out : ../src/xbsstep.f $(ILIB)
	$(OUTRULE)

xcaldat.out : ../src/xcaldat.f $(ILIB)
	$(OUTRULE)

xchder.out : ../src/xchder.f $(ILIB)
	$(OUTRULE)

xchebev.out : ../src/xchebev.f $(ILIB)
	$(OUTRULE)

xchebft.out : ../src/xchebft.f $(ILIB)
	$(OUTRULE)

xchebpc.out : ../src/xchebpc.f $(ILIB)
	$(OUTRULE)

xchint.out : ../src/xchint.f $(ILIB)
	$(OUTRULE)

xcholsl.out : ../src/xcholsl.f $(ILIB)
	$(OUTRULE)

xchsone.out : ../src/xchsone.f $(ILIB)
	$(OUTRULE)

xchstwo.out : ../src/xchstwo.f $(ILIB)
	$(OUTRULE)

xcisi.out : ../src/xcisi.f $(ILIB)
	$(OUTRULE)

xcntab1.out : ../src/xcntab1.f $(ILIB)
	$(OUTRULE)

xcntab2.out : ../src/xcntab2.f $(ILIB)
	$(OUTRULE)

xconvlv.out : ../src/xconvlv.f $(ILIB)
	$(OUTRULE)

xcorrel.out : ../src/xcorrel.f $(ILIB)
	$(OUTRULE)

xcosft1.out : ../src/xcosft1.f $(ILIB)
	$(OUTRULE)

xcosft2.out : ../src/xcosft2.f $(ILIB)
	$(OUTRULE)

xcovsrt.out : ../src/xcovsrt.f $(ILIB)
	$(OUTRULE)

xcrank.out : ../src/xcrank.f $(ILIB)
	$(OUTRULE)

xcyclic.out : ../src/xcyclic.f $(ILIB)
	$(OUTRULE)

xdawson.out : ../src/xdawson.f $(ILIB)
	$(OUTRULE)

xdbrent.out : ../src/xdbrent.f $(ILIB)
	$(OUTRULE)

xddpoly.out : ../src/xddpoly.f $(ILIB)
	$(OUTRULE)

xdecchk.out : ../src/xdecchk.f $(ILIB)
	$(OUTRULE)

xdf1dim.out : ../src/xdf1dim.f $(ILIB)
	$(OUTRULE)

xdfpmin.out : ../src/xdfpmin.f $(ILIB)
	$(OUTRULE)

xdfridr.out : ../src/xdfridr.f $(ILIB)
	$(OUTRULE)

xdftint.out : ../src/xdftint.f $(ILIB)
	$(OUTRULE)

xeclass.out : ../src/xeclass.f $(ILIB)
	$(OUTRULE)

xeclazz.out : ../src/xeclazz.f $(ILIB)
	$(OUTRULE)

xei.out : ../src/xei.f $(ILIB)
	$(OUTRULE)

xeigsrt.out : ../src/xeigsrt.f $(ILIB)
	$(OUTRULE)

xelle.out : ../src/xelle.f $(ILIB)
	$(OUTRULE)

xellf.out : ../src/xellf.f $(ILIB)
	$(OUTRULE)

xellpi.out : ../src/xellpi.f $(ILIB)
	$(OUTRULE)

xelmhes.out : ../src/xelmhes.f $(ILIB)
	$(OUTRULE)

xerf.out : ../src/xerf.f $(ILIB)
	$(OUTRULE)

xerfc.out : ../src/xerfc.f $(ILIB)
	$(OUTRULE)

xerfcc.out : ../src/xerfcc.f $(ILIB)
	$(OUTRULE)

xeulsum.out : ../src/xeulsum.f $(ILIB)
	$(OUTRULE)

xevlmem.out : ../src/xevlmem.f $(ILIB)
	$(OUTRULE)

xexpdev.out : ../src/xexpdev.f $(ILIB)
	$(OUTRULE)

xexpint.out : ../src/xexpint.f $(ILIB)
	$(OUTRULE)

xf1dim.out : ../src/xf1dim.f $(ILIB)
	$(OUTRULE)

xfactln.out : ../src/xfactln.f $(ILIB)
	$(OUTRULE)

xfactrl.out : ../src/xfactrl.f $(ILIB)
	$(OUTRULE)

xfasper.out : ../src/xfasper.f $(ILIB)
	$(OUTRULE)

xfgauss.out : ../src/xfgauss.f $(ILIB)
	$(OUTRULE)

xfit.out : ../src/xfit.f $(ILIB)
	$(OUTRULE)

xfitexy.out : ../src/xfitexy.f $(ILIB)
	$(OUTRULE)

xfixrts.out : ../src/xfixrts.f $(ILIB)
	$(OUTRULE)

xfleg.out : ../src/xfleg.f $(ILIB)
	$(OUTRULE)

xflmoon.out : ../src/xflmoon.f $(ILIB)
	$(OUTRULE)

xfour1.out : ../src/xfour1.f $(ILIB)
	$(OUTRULE)

xfourfs.out : ../src/xfourfs.f $(ILIB)
	$(OUTRULE)

xfourn.out : ../src/xfourn.f $(ILIB)
	$(OUTRULE)

xfpoly.out : ../src/xfpoly.f $(ILIB)
	$(OUTRULE)

xfred2.out : ../src/xfred2.f $(ILIB)
	$(OUTRULE)

xfredin.out : ../src/xfredin.f $(ILIB)
	$(OUTRULE)

xfrenel.out : ../src/xfrenel.f $(ILIB)
	$(OUTRULE)

xfrprmn.out : ../src/xfrprmn.f $(ILIB)
	$(OUTRULE)

xftest.out : ../src/xftest.f $(ILIB)
	$(OUTRULE)

xgamdev.out : ../src/xgamdev.f $(ILIB)
	$(OUTRULE)

xgammln.out : ../src/xgammln.f $(ILIB)
	$(OUTRULE)

xgammp.out : ../src/xgammp.f $(ILIB)
	$(OUTRULE)

xgammq.out : ../src/xgammq.f $(ILIB)
	$(OUTRULE)

xgasdev.out : ../src/xgasdev.f $(ILIB)
	$(OUTRULE)

xgaucof.out : ../src/xgaucof.f $(ILIB)
	$(OUTRULE)

xgauher.out : ../src/xgauher.f $(ILIB)
	$(OUTRULE)

xgaujac.out : ../src/xgaujac.f $(ILIB)
	$(OUTRULE)

xgaulag.out : ../src/xgaulag.f $(ILIB)
	$(OUTRULE)

xgauleg.out : ../src/xgauleg.f $(ILIB)
	$(OUTRULE)

xgaussj.out : ../src/xgaussj.f $(ILIB)
	$(OUTRULE)

xgcf.out : ../src/xgcf.f $(ILIB)
	$(OUTRULE)

xgolden.out : ../src/xgolden.f $(ILIB)
	$(OUTRULE)

xgser.out : ../src/xgser.f $(ILIB)
	$(OUTRULE)

xhpsel.out : ../src/xhpsel.f $(ILIB)
	$(OUTRULE)

xhpsort.out : ../src/xhpsort.f $(ILIB)
	$(OUTRULE)

xhqr.out : ../src/xhqr.f $(ILIB)
	$(OUTRULE)

xhuffman.out : ../src/xhuffman.f $(ILIB)
	$(OUTRULE)

xhunt.out : ../src/xhunt.f $(ILIB)
	$(OUTRULE)

xhypgeo.out : ../src/xhypgeo.f $(ILIB)
	$(OUTRULE)

xicrc.out : ../src/xicrc.f $(ILIB)
	$(OUTRULE)

xigray.out : ../src/xigray.f $(ILIB)
	$(OUTRULE)

xindexx.out : ../src/xindexx.f $(ILIB)
	$(OUTRULE)

xirbit1.out : ../src/xirbit1.f $(ILIB)
	$(OUTRULE)

xirbit2.out : ../src/xirbit2.f $(ILIB)
	$(OUTRULE)

xjacobi.out : ../src/xjacobi.f $(ILIB)
	$(OUTRULE)

xjulday.out : ../src/xjulday.f $(ILIB)
	$(OUTRULE)

xkendl1.out : ../src/xkendl1.f $(ILIB)
	$(OUTRULE)

xkendl2.out : ../src/xkendl2.f $(ILIB)
	$(OUTRULE)

xks2d1s.out : ../src/xks2d1s.f $(ILIB)
	$(OUTRULE)

xks2d2s.out : ../src/xks2d2s.f $(ILIB)
	$(OUTRULE)

xksone.out : ../src/xksone.f $(ILIB)
	$(OUTRULE)

xkstwo.out : ../src/xkstwo.f $(ILIB)
	$(OUTRULE)

xlaguer.out : ../src/xlaguer.f $(ILIB)
	$(OUTRULE)

xlfit.out : ../src/xlfit.f $(ILIB)
	$(OUTRULE)

xlinbcg.out : ../src/xlinbcg.f $(ILIB)
	$(OUTRULE)

xlinmin.out : ../src/xlinmin.f $(ILIB)
	$(OUTRULE)

xlocate.out : ../src/xlocate.f $(ILIB)
	$(OUTRULE)

xlubksb.out : ../src/xlubksb.f $(ILIB)
	$(OUTRULE)

xludcmp.out : ../src/xludcmp.f $(ILIB)
	$(OUTRULE)

xmachar.out : ../src/xmachar.f $(ILIB)
	$(OUTRULE)

xmedfit.out : ../src/xmedfit.f $(ILIB)
	$(OUTRULE)

xmemcof.out : ../src/xmemcof.f $(ILIB)
	$(OUTRULE)

xmgfas.out : ../src/xmgfas.f $(ILIB)
	$(OUTRULE)

xmglin.out : ../src/xmglin.f $(ILIB)
	$(OUTRULE)

xmidpnt.out : ../src/xmidpnt.f $(ILIB)
	$(OUTRULE)

xmiser.out : ../src/xmiser.f $(ILIB)
	$(OUTRULE)

xmmid.out : ../src/xmmid.f $(ILIB)
	$(OUTRULE)

xmnbrak.out : ../src/xmnbrak.f $(ILIB)
	$(OUTRULE)

xmnewt.out : ../src/xmnewt.f $(ILIB)
	$(OUTRULE)

xmoment.out : ../src/xmoment.f $(ILIB)
	$(OUTRULE)

xmppi.out : ../src/xmppi.f $(ILIB)
	$(OUTRULE)

xmprove.out : ../src/xmprove.f $(ILIB)
	$(OUTRULE)

xmrqcof.out : ../src/xmrqcof.f $(ILIB)
	$(OUTRULE)

xmrqmin.out : ../src/xmrqmin.f $(ILIB)
	$(OUTRULE)

xnewt.out : ../src/xnewt.f $(ILIB)
	$(OUTRULE)

xodeint.out : ../src/xodeint.f $(ILIB)
	$(OUTRULE)

xorthog.out : ../src/xorthog.f $(ILIB)
	$(OUTRULE)

xpade.out : ../src/xpade.f $(ILIB)
	$(OUTRULE)

xpccheb.out : ../src/xpccheb.f $(ILIB)
	$(OUTRULE)

xpcshft.out : ../src/xpcshft.f $(ILIB)
	$(OUTRULE)

xpearsn.out : ../src/xpearsn.f $(ILIB)
	$(OUTRULE)

xperiod.out : ../src/xperiod.f $(ILIB)
	$(OUTRULE)

xpiksr2.out : ../src/xpiksr2.f $(ILIB)
	$(OUTRULE)

xpiksrt.out : ../src/xpiksrt.f $(ILIB)
	$(OUTRULE)

xplgndr.out : ../src/xplgndr.f $(ILIB)
	$(OUTRULE)

xpoidev.out : ../src/xpoidev.f $(ILIB)
	$(OUTRULE)

xpolcoe.out : ../src/xpolcoe.f $(ILIB)
	$(OUTRULE)

xpolcof.out : ../src/xpolcof.f $(ILIB)
	$(OUTRULE)

xpoldiv.out : ../src/xpoldiv.f $(ILIB)
	$(OUTRULE)

xpolin2.out : ../src/xpolin2.f $(ILIB)
	$(OUTRULE)

xpolint.out : ../src/xpolint.f $(ILIB)
	$(OUTRULE)

xpowell.out : ../src/xpowell.f $(ILIB)
	$(OUTRULE)

xpredic.out : ../src/xpredic.f $(ILIB)
	$(OUTRULE)

xprobks.out : ../src/xprobks.f $(ILIB)
	$(OUTRULE)

xpsdes.out : ../src/xpsdes.f $(ILIB)
	$(OUTRULE)

xpzextr.out : ../src/xpzextr.f $(ILIB)
	$(OUTRULE)

xqgaus.out : ../src/xqgaus.f $(ILIB)
	$(OUTRULE)

xqrdcmp.out : ../src/xqrdcmp.f $(ILIB)
	$(OUTRULE)

xqromb.out : ../src/xqromb.f $(ILIB)
	$(OUTRULE)

xqromo.out : ../src/xqromo.f $(ILIB)
	$(OUTRULE)

xqroot.out : ../src/xqroot.f $(ILIB)
	$(OUTRULE)

xqrsolv.out : ../src/xqrsolv.f $(ILIB)
	$(OUTRULE)

xqrupdt.out : ../src/xqrupdt.f $(ILIB)
	$(OUTRULE)

xqsimp.out : ../src/xqsimp.f $(ILIB)
	$(OUTRULE)

xqtrap.out : ../src/xqtrap.f $(ILIB)
	$(OUTRULE)

xquad3d.out : ../src/xquad3d.f $(ILIB)
	$(OUTRULE)

xran.out : ../src/xran.f $(ILIB)
	$(OUTRULE)

xran4.out : ../src/xran4.f $(ILIB)
	$(OUTRULE)

xrank.out : ../src/xrank.f $(ILIB)
	$(OUTRULE)

xratint.out : ../src/xratint.f $(ILIB)
	$(OUTRULE)

xratlsq.out : ../src/xratlsq.f $(ILIB)
	$(OUTRULE)

xrc.out : ../src/xrc.f $(ILIB)
	$(OUTRULE)

xrd.out : ../src/xrd.f $(ILIB)
	$(OUTRULE)

xrealft.out : ../src/xrealft.f $(ILIB)
	$(OUTRULE)

xrf.out : ../src/xrf.f $(ILIB)
	$(OUTRULE)

xrj.out : ../src/xrj.f $(ILIB)
	$(OUTRULE)

xrk4.out : ../src/xrk4.f $(ILIB)
	$(OUTRULE)

xrkdumb.out : ../src/xrkdumb.f $(ILIB)
	$(OUTRULE)

xrkqs.out : ../src/xrkqs.f $(ILIB)
	$(OUTRULE)

xrlft3.out : ../src/xrlft3.f $(ILIB)
	$(OUTRULE)

xrofunc.out : ../src/xrofunc.f $(ILIB)
	$(OUTRULE)

xrtbis.out : ../src/xrtbis.f $(ILIB)
	$(OUTRULE)

xrtflsp.out : ../src/xrtflsp.f $(ILIB)
	$(OUTRULE)

xrtnewt.out : ../src/xrtnewt.f $(ILIB)
	$(OUTRULE)

xrtsafe.out : ../src/xrtsafe.f $(ILIB)
	$(OUTRULE)

xrtsec.out : ../src/xrtsec.f $(ILIB)
	$(OUTRULE)

xrzextr.out : ../src/xrzextr.f $(ILIB)
	$(OUTRULE)

xsavgol.out : ../src/xsavgol.f $(ILIB)
	$(OUTRULE)

xscrsho.out : ../src/xscrsho.f $(ILIB)
	$(OUTRULE)

xselect.out : ../src/xselect.f $(ILIB)
	$(OUTRULE)

xselip.out : ../src/xselip.f $(ILIB)
	$(OUTRULE)

xshell.out : ../src/xshell.f $(ILIB)
	$(OUTRULE)

xsimplx.out : ../src/xsimplx.f $(ILIB)
	$(OUTRULE)

xsimpr.out : ../src/xsimpr.f $(ILIB)
	$(OUTRULE)

xsinft.out : ../src/xsinft.f $(ILIB)
	$(OUTRULE)

xsncndn.out : ../src/xsncndn.f $(ILIB)
	$(OUTRULE)

xsobseq.out : ../src/xsobseq.f $(ILIB)
	$(OUTRULE)

xsor.out : ../src/xsor.f $(ILIB)
	$(OUTRULE)

xsort.out : ../src/xsort.f $(ILIB)
	$(OUTRULE)

xsort2.out : ../src/xsort2.f $(ILIB)
	$(OUTRULE)

xsort3.out : ../src/xsort3.f $(ILIB)
	$(OUTRULE)

xspctrm.out : ../src/xspctrm.f $(ILIB)
	$(OUTRULE)

xspear.out : ../src/xspear.f $(ILIB)
	$(OUTRULE)

xsphbes.out : ../src/xsphbes.f $(ILIB)
	$(OUTRULE)

xsplie2.out : ../src/xsplie2.f $(ILIB)
	$(OUTRULE)

xsplin2.out : ../src/xsplin2.f $(ILIB)
	$(OUTRULE)

xspline.out : ../src/xspline.f $(ILIB)
	$(OUTRULE)

xsplint.out : ../src/xsplint.f $(ILIB)
	$(OUTRULE)

xsprsax.out : ../src/xsprsax.f $(ILIB)
	$(OUTRULE)

xsprsin.out : ../src/xsprsin.f $(ILIB)
	$(OUTRULE)

xsprspm.out : ../src/xsprspm.f $(ILIB)
	$(OUTRULE)

xsprstm.out : ../src/xsprstm.f $(ILIB)
	$(OUTRULE)

xsprstp.out : ../src/xsprstp.f $(ILIB)
	$(OUTRULE)

xsprstx.out : ../src/xsprstx.f $(ILIB)
	$(OUTRULE)

xstifbs.out : ../src/xstifbs.f $(ILIB)
	$(OUTRULE)

xstiff.out : ../src/xstiff.f $(ILIB)
	$(OUTRULE)

xstoerm.out : ../src/xstoerm.f $(ILIB)
	$(OUTRULE)

xsvbksb.out : ../src/xsvbksb.f $(ILIB)
	$(OUTRULE)

xsvdfit.out : ../src/xsvdfit.f $(ILIB)
	$(OUTRULE)

xsvdvar.out : ../src/xsvdvar.f $(ILIB)
	$(OUTRULE)

xtoeplz.out : ../src/xtoeplz.f $(ILIB)
	$(OUTRULE)

xtptest.out : ../src/xtptest.f $(ILIB)
	$(OUTRULE)

xtqli.out : ../src/xtqli.f $(ILIB)
	$(OUTRULE)

xtrapzd.out : ../src/xtrapzd.f $(ILIB)
	$(OUTRULE)

xtred2.out : ../src/xtred2.f $(ILIB)
	$(OUTRULE)

xtridag.out : ../src/xtridag.f $(ILIB)
	$(OUTRULE)

xttest.out : ../src/xttest.f $(ILIB)
	$(OUTRULE)

xtutest.out : ../src/xtutest.f $(ILIB)
	$(OUTRULE)

xtwofft.out : ../src/xtwofft.f $(ILIB)
	$(OUTRULE)

xvander.out : ../src/xvander.f $(ILIB)
	$(OUTRULE)

xvegas.out : ../src/xvegas.f $(ILIB)
	$(OUTRULE)

xvoltra.out : ../src/xvoltra.f $(ILIB)
	$(OUTRULE)

xwt1.out : ../src/xwt1.f $(ILIB)
	$(OUTRULE)

xwtn.out : ../src/xwtn.f $(ILIB)
	$(OUTRULE)

xzbrac.out : ../src/xzbrac.f $(ILIB)
	$(OUTRULE)

xzbrak.out : ../src/xzbrak.f $(ILIB)
	$(OUTRULE)

xzbrent.out : ../src/xzbrent.f $(ILIB)
	$(OUTRULE)

xzrhqr.out : ../src/xzrhqr.f $(ILIB)
	$(OUTRULE)

xzriddr.out : ../src/xzriddr.f $(ILIB)
	$(OUTRULE)

xzroots.out : ../src/xzroots.f $(ILIB)
	$(OUTRULE)
