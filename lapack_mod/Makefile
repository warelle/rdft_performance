all:
	gfortran -Wall -O3 -c dgetrfw.f dgetrfw2.f zgetrfw.f zgetrfw2.f dlamch.f ieeeck.f ilaenv.f iparmq.f lsame.f xerbla.f
	ar rc mod_lu.a dgetrfw.o dgetrfw2.o zgetrfw.o zgetrfw2.o dlamch.o ieeeck.o ilaenv.o iparmq.o lsame.o

clean:
	rm *.o *~
