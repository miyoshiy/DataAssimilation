C     Functions which implement some bit functions in terms of
C     other intrinsics
C     This file is provided for convenience. It is simply one way of dealing
C     with bit function name differences.

      LOGICAL function btest(iword,ibitnum)
      INTEGER iword,ibitnum,and,rshift
      if(and(rshift(iword,ibitnum),1).eq.1) then
         btest=.true.
      else
         btest=.false.
      endif
      return
      END

      INTEGER function ibset(iword,ibitnum)
      INTEGER iword,ibitnum,or,lshift
      ibset=or(iword,lshift(1,ibitnum))
      return
      END

      INTEGER function ibclr(iword,ibitnum)
      INTEGER iword,ibitnum,and,not,lshift
      ibclr=and(iword,not(lshift(1,ibitnum)))
      return
      END

      INTEGER function ieor(iword1,iword2)
      INTEGER iword1,iword2,xor
      ieor=xor(iword1,iword2)
      return
      END

      INTEGER function ior(iword1,iword2)
      INTEGER iword1,iword2,or
      ior=or(iword1,iword2)
      return
      END

      INTEGER function iand(iword1,iword2)
      INTEGER iword1,iword2,and
      iand=and(iword1,iword2)
      return
      END

      INTEGER function ishft(iword,nbits)
      INTEGER nbits,iword,lshift,rshift
      if (nbits.ge.0) then
        ishft=lshift(iword,nbits)
      else
        ishft=rshift(iword,-nbits)
      endif
      return
      END

C     Function which fixes a problem in ichar in certain machines.
C     This returns the range [0, 255], even though ichar
C     may return from range [-128, 127].
      INTEGER function iichar(ch)
      INTEGER ichar
      CHARACTER ch
      iichar=ichar(ch)
      if (iichar.lt.0) iichar=iichar+256
      return
      END
