#!/bin/bash
#=====================================================
#Fontes
bold=$(tput bold)
normal=$(tput sgr0)
red=$(tput setaf 1)
green=$(tput setaf 2)
yellow=$(tput setaf 3)
blue=$(tput setaf 4)
pink=$(tput setaf 5)
blue2=$(tput setaf 6)
grey=$(tput setaf 7)
white=$(tput setaf 8)
reset=$(tput sgr0)
#=====================================================
#CABEÇALHO DO PROGRAMA
echo "${blue}${bold}==========================================================================================${normal}"
echo "       ${red}UTFPR - UNIVERSIDADE TECNOLÓGICA FEDERAL DO PARANÁ - DEPARTAMENTO DE FÍSICA ${normal}"
echo "   ${red}${bold}CONVERSOR DE RESULTADOS DO P.R.A.I.A. (IRAF) PARA O A.T.O.S. (AUDELA) - Versão 0.4b${normal}"
echo ""
echo "Desenvolvido por ${blue}${bold}RICARDO CEZAR VOLERT${normal} em 01-oct-2015"
echo "E-mail: ${bold}ricardovolert@alunos.utfpr.edu.br${normal}"
echo "Atualizado em 13-nov-2015"
echo "(Não funciona em Ubuntu 17.04 - Verificado em 28-out-2017)"
echo "     ${green}Agradecimentos pela ajuda no desenvolvimento à: Felipe Braga-Ribas e Douglas Mancini${normal}"
echo ""
echo "    ${blue}${bold}+------------------------------------------------------------------------------------+${normal}"
echo "    ${blue}${bold}|${normal}                                                                                    ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal}                                ${blue}${bold}RESUMO DO PROGRAMA:${normal}                                 ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal}                                                                                    ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal} O programa filtra colunas dos arquivos de fotometria e curva de luz gerado pelo    ${blue}${bold}|${normal}" 
echo "    ${blue}${bold}|${normal} ${bold}P.R.A.I.A.${normal} (Software ${bold}${pink}IRAF${normal}) e gera novos arquivos, compatíveis com o ${bold}A.T.O.S.${normal}       ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal} (Software ${bold}${pink}AUDELA${normal}) e gera um gráfico em formato ${bold}PNG${normal} para visualização prévia da     ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal} curva a partir do ${bold}GNUPLOT${normal} que deverá estar instalado.                              ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal}                                                                                    ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal} ${bold}AVISOS${normal}:                                                                            ${blue}${bold}|${normal}" 
echo "    ${blue}${bold}|${normal}     1) Este programa funciona juntamente com uso do ${bold}GFORTRAN${normal}, ${bold}GNUPLOT${normal} e ${bold}SHOTWELL${normal}.  ${blue}${bold}|${normal}"
echo "    ${blue}${bold}|${normal}     2) Desenvolvido e Funcionando em distribuição ${bold}LINUX OpenSUSE 13.1 - 64 bits${normal}.   ${blue}${bold}|${normal}"
echo "    ${blue}${bold}+------------------------------------------------------------------------------------+${normal}"
echo ""
#=====================================================


curve=$(ls -c curve_*.dat | wc -l)
photometry=$(ls -c photometry_*.dat | wc -l)
unicofile=1
opcaodetrabalho=0

if (($photometry == "1" )) || (($photometry == "2" )); then
	photometry=1
else
	if (($photometry != "1" )) || (($photometry != "2" )); then
		photometry=0
	fi
fi
if (($curve == "1" )) || (($curve == "2" )); then
	curve=1
else
	if (($curve != "1" )) || (($curve != "2" )); then
		curve=0
	fi
fi


if (("$curve" >= "$unicofile")) && (("$photometry" >= "$unicofile"))
then
	echo ""
	echo "Você possuí tanto arquivos de fotometria quanto de curvas de luz, escolha qual opção deseja:"
	echo ""
	echo "            ${bold}(1)${normal} Curva de Luz (curve_*.dat)  -  ${bold}(2)${normal} Fotometria (photometry_*.dat)"
	echo ""
	while (("$opcaodetrabalho" > "2")) || (("$opcaodetrabalho" < "1")); do
		echo "Por favor insira o número apenas de sua opção escolhida, ou 'q', 'Q', 'exit' para sair:"
		read opcaodetrabalho
		if (("$opcaodetrabalho" == "q")) || (("$opcaodetrabalho" == "Q")) || (("$opcaodetrabalho" == "exit"))
		then
			exit
		fi
	done
	#echo $opcaodetrabalho
else
	if (("$curve" >= "$unicofile"))
	then
		opcaodetrabalho=1
		#echo $opcaodetrabalho
	else
		echo "${red}${bold}"
		echo "   * ALERTA: Não foi detectado arquivo curve_*.dat${normal}"
	fi
	if (("$photometry" >= "$unicofile"))
	then
		opcaodetrabalho=2
		#echo $opcaodetrabalho
	else
		echo "${red}${bold}"
		echo "   * ALERTA: Não foi detectado arquivo photometry_*.dat${normal}"
	fi
	if (("$curve" == "0")) && (("$photometry" == "0"))
	then
		echo "${red}${bold}"
		echo "   * AVISO: Como o sistema não encontrou nenhum dos formatos acima relatado, o programa" 
		echo "     ${blue}IRAF2Audela${red} está sendo encerrado${normal}"
		exit
	fi
fi
#=====================================================================
#Se foi escolhido a opção curva de luz (1)
if (("$opcaodetrabalho" == "1"))
then
#=====================================================================
#Procura quantas linhas dará
#gawk '{ sum += $1 }; END { print sum }' curve_*.dat > contalinhas.txt
echo "Por Favor, insira o numero de linhas do arquivo curve:"
read totaldelinhascurve
echo "$totaldelinhascurve" > contalinhas.txt

#Cria Fit File
gawk '{ sum2 += $1 }; END { print (sum2 / 5)}' curve_*.dat > intervalo_de_tempo_gnuplot.txt
contouaslinhas=$(<intervalo_de_tempo_gnuplot.txt)
echo "1" > tics_number_gnuplot.txt
echo "$contouaslinhas * 1" | bc -l >> tics_number_gnuplot.txt
echo "$contouaslinhas * 2" | bc -l >> tics_number_gnuplot.txt
echo "$contouaslinhas * 3" | bc -l >> tics_number_gnuplot.txt
echo "$contouaslinhas * 4" | bc -l >> tics_number_gnuplot.txt
gawk '{ sum += $1 }; END { print sum }' curve_*.dat >> tics_number_gnuplot.txt

#Arredonda os valores
cat <<EOF >arrendadomento_de_quadros.f90
Program Arredondamento
	Implicit None
! ======================================================================================
	Integer :: contaarredondador	
	Real(Kind = 8), dimension(6) :: arredondar
	OPEN(Unit=21, file="tics_number_gnuplot.txt", status="Old")
	OPEN(Unit=22, file="tics_number_gnuplot_corrigido.txt", status="Replace")
	do contaarredondador=1, 6
		read(21,*) arredondar(contaarredondador)
		Write(22,*) ceiling(arredondar(contaarredondador))		 
	end do
	close (Unit=22)	
	close (Unit=21)
End Program Arredondamento
EOF

#Encerra Fit File
	gfortran 'arrendadomento_de_quadros.f90' -o 'arrendadomento_de_quadros'
	./arrendadomento_de_quadros


for filename in $(ls -c curve_*.dat | head -1);
do
#=============================================================
#Insere o seguinte código-fonte em um arquivo do Fortran Gerado como Padrão
cat <<EOF >PRAIA-results_to_AudeLa.f90
Program PRAIAtoAudeLA
	Implicit None
! ======================================================================================

	Real(Kind = 8) :: U02, U03, U04, U05, U06, U07, U08, U09, U10, U11, U13, b
	DOUBLE PRECISION, dimension(100000) :: JDay
	Real(Kind = 8), dimension(100000) :: BrilhoObjeto , BrilhoRefI 
	Integer :: i, j, a, contalinhas
	CHARACTER(20) :: NomeDoFit
	CHARACTER(13) :: TimeISOa
	CHARACTER(2) :: TimeISOb
	CHARACTER(7) :: TimeISOc
	CHARACTER(27) :: TimeISO
	 
! ======================================================================================
	Write(*,*) " "
	Write(*,*) " "

! ======================================================================================
!Arquivo: arquivo_velho
TimeISO = ''
! Leitura nos arquivos de dados
	OPEN(Unit=15, file="contalinhas.txt", status="Old")
		read(15,*) contalinhas
		!contalinhas = 251 
		!write(*,*) "Linhas de Texto", contalinhas
	close (Unit=15)

	OPEN(Unit=2, file="curve_OBJETOESTUDADO.dat", status="Old")
	! Gerando arquivo corrigido
	  OPEN(Unit=3, file="RESULT_OBJETOESTUDADO_.00000.csv", status="Replace")
	  OPEN(Unit=99, file="RESULT_OBJETOESTUDADO_.time", status="Replace")
	  OPEN(Unit=98, file="RESULT_.freemat.00000.dat", status="Replace")
	  OPEN(Unit=97, file="RESULT_deltatime.dat", status="Replace")
		Write(3,*) "# ** atos - Audela - Linux  * "
		Write(99,*) "# ** atos - Audela - Linux  * "
		Write(3,*) "# FPS 25"
		Write(99,*) "# FPS 25"		
		Write(3,*) "idframe,jd,dateiso,obj_fint     ,obj_pixmax   ,obj_intensite, &
obj_sigmafond,obj_snint    ,obj_snpx     ,obj_delta    ,obj_xpos,obj_ypos,obj_xfwhm,obj_yfwhm,&
ref_fint     ,ref_pixmax   ,ref_intensite,ref_sigmafond,ref_snint    ,ref_snpx     ,ref_delta    ,&
ref_xpos,ref_ypos,ref_xfwhm,ref_yfwhm,img_intmin ,img_intmax,img_intmoy,img_sigma ,img_xsize,img_ysize"
		Write(99,*) "idframe, jd, dateiso, verif, ocr, interpol"
		do i=1, contalinhas

			! a, b, JDay(i), U02, U04, U05, U06, U07, U08, U09, U10, BrilhoObjeto(i), BrilhoRefI(i)

	  		read(2,*) a, b, JDay(i), U02, U04, U05, U06, U07, U08, U09, &
U10, BrilhoObjeto(i), BrilhoRefI(i)
!	  		write(*,*) a, b, JDay(i), U02, U04, U05, U06, U07, U08, U09, &
!U10, BrilhoObjeto(i), BrilhoRefI(i), " "
			Call jd2greg(JDay(i))
			
			OPEN(Unit=11, file="auxiliar2.dat", status="Old")
				do j=1,1
					read(11,*) TimeISOa, TimeISOb, TimeISOc
					!concatenação
					TimeISO = TimeISOa //':'// TimeISOb //':'// TimeISOc
				end do
			close (Unit=11)

			write(97,*) TimeISO(12:13), "    " , TimeISO(15:16), "    " , TimeISO(18:24)

			write(3,10010) i, JDay(i), TimeISO, BrilhoObjeto(i), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, BrilhoRefI(i), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			write(99,10030) i, JDay(i), TimeISO, 1, 1, 1

			write(98,10020) i, BrilhoObjeto(i), BrilhoRefI(i)
	   	end do
	close(97)	
	close(98)	
	close(99)	  
	close(3)
	close(Unit = 2)



	!Depois de tudo feito
	OPEN(Unit=11, file="auxiliar.dat", status="Old")
	close (Unit=11, status='delete')
	OPEN(Unit=11, file="auxiliar2.dat", status="Old")
	close (Unit=11, status='delete')

	Write(*,*) " "
	Write(*,*) "Resultados de Fotometria foram salvos com sucesso!"
	Write(*,*) "=========================================================================================="

	Write(*,*) " "
	Write(*,*) " "

! ======================================================================================
! Formatos de Saida
10010 FORMAT(I5.1, ',', F16.8, ',', A25, ',', F15.4, ',', I3.1, &
',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', F15.4, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1)

10020 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4)

10030 FORMAT(I5.1, ',', F16.8, ',', A25, ',', I3.1, ',', I3.1, ',', I3.1)

10000 FORMAT(' ', I5.1, ',', F16.8, ',', A25, ',', F10.4, ',', I3.1, ',', A30) !Relatorio Final
10001 FORMAT('i = ',I3.1) !para inteiro
10002 FORMAT(5x,'JD(i) = ',F16.8) !para float de JD
10003 FORMAT(5x,'BrilhoObjeto(i) = ',F9.4) !float de brilho
10004 FORMAT(5x,'Nome do Fit ='A30) !texto

! ======================================================================================

End Program PRAIAtoAudeLA
! ======================================================================================
! Subrotinas que convertem uma data JD para ISO (Julian Day to a Gregorian Date)

SUBROUTINE jd2greg (julianTime)
	!DOUBLE PRECISION :: julianTime
	!julianTime = 2*julianTime
	!write(*,20002) julianTime
	!20002 FORMAT(5x,'JDay(i) = ',F16.8) !para float de JD

      IMPLICIT REAL *8 (A-H,O-Z)
      DOUBLE PRECISION :: julianTime
      character*24 ler, julianTimeText
      hmsgms(ai,aj,aa)=ai+aj/60.d0+aa/3600.d0
      dj2000=2451544.5d0
      dj1=2400000.5D0
      idim=250

      julianTimeText = ''
    open(4,file="auxiliar.dat")
       write(4,*) julianTime
    close(4)

    open(4,file="auxiliar.dat")
       read(4, *) julianTimeText
       !write(*,*) A_char
    close(4)

      ler=julianTimeText !2456977.23464133'
!'
      read (ler,*) data
      iflag=1
      if (data.lt.2500.d0) iflag=2
      if (data.gt.25.d0.and.data.lt.2500.d0) iflag=3
      if (iflag.ne.1) go to 30
 15   continue
      djm=data-dj1
      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,jjj)
      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0
      frano=(data-dj2000)/365.25d0+2000.d0
      ler=''
!      write (ler,20) iuth,iutm,sut,iutdia,iutmes,iutano,data,frano
! 20   format(2(i2.2,':'),f7.4,':',2(i2.2,':'),i4.4,':',f18.10,':',
!     ?f15.10)

!      write (ler,20) iuth,iutm,sut,iutdia,iutmes,iutano
    open(666,file="auxiliar2.dat")
 !      write (666,30001) iutano,iutmes,iutdia,iuth,iutm,sut
 !30001   format(i4.4,'-',1(i2.2,'-'),i2.2,'T', 2(i2.2,':'),f7.4)
   ! close(666)


      write (ler,20) iutano,iutmes,iutdia,iuth,iutm,sut
!	write (*,20) iutano,iutmes,iutdia,iuth,iutm,sut
 20   format(i4.4,'-',1(i2.2,'-'),i2.2,'T', 2(i2.2,':'),f7.4)
!20   format(i4.4,'-',1(i2.2,'-'),i2.2,'T', 2(i2.2,':'),f7.4)


      do i=idim,1,-1
      if (ler(i:i).ne.' ') go to 25
      enddo
 25   ii=i

      do i=1,ii-1
      if (ler(i:i).eq.' ') ler(i:i)='0'
      enddo

      do i=1,ii-1
      if (ler(i:i).eq.':') ler(i:i)=' '
      enddo

!      write (*,*) ler
	write (666,*) ler
	close(666)
      go to 100

 30   continue

      if (iflag.eq.3) then
      data=365.25d0*(data-2000.d0)+dj2000
      go to 15
      endif

      read (ler,*) uth,utm,sut,utdia,iutmes,iutano

      fd=hmsgms(uth,utm,sut)/24.d0

      iutdia=utdia

      fdd=utdia-iutdia

      fd=fd+fdd

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,jjj)


      djm=djm+fd
            
      data=djm+djm0


      go to 15

 100  continue




END SUBROUTINE jd2greg



      SUBROUTINE iau_CAL2JD ( IY, IM, ID, DJM0, DJM, J )

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION DJM0, DJM
      INTEGER J, MY, IYPMY

!  Earliest year allowed (4800BC)
      INTEGER IYMIN
      PARAMETER ( IYMIN = -4799 )

!  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Preset status.
      J = 0

!  Validate year.
      IF ( IY.LT.IYMIN ) THEN
         J = -1
      ELSE

!     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN

!        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 ) MTAB(2) = 28

!        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J = -3

!        Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( (1461*(IYPMY+4800))/4 + (367*(IM-2-(12*MY)))/12 - (3*((IYPMY+4900)/100))/4 + ID-2432076)

!        Bad month
         ELSE
            J = -2
         END IF
      END IF


      END


      SUBROUTINE iau_jd2cal ( DJ1, DJ2, IY, IM, ID, FD, J )

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2
      INTEGER IY, IM, ID
      DOUBLE PRECISION FD
      INTEGER J

!  Minimum and maximum allowed JD
      DOUBLE PRECISION DJMIN, DJMAX
      PARAMETER ( DJMIN = -68569.5D0, DJMAX = 1D9 )

      INTEGER JD, L, N, I
      DOUBLE PRECISION DJ, D1, D2, F1, F2, F, D

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Check if date is acceptable.
      DJ = DJ1 + DJ2
      IF ( DJ.LT.DJMIN .OR. DJ.GT.DJMAX ) THEN
         J = -1
      ELSE
         J = 0

!     Copy the date, big then small, and re-align to midnight.
         IF ( DJ1 .GE. DJ2 ) THEN
            D1 = DJ1
            D2 = DJ2
         ELSE
            D1 = DJ2
            D2 = DJ1
         END IF
         D2 = D2 - 0.5D0

!     Separate day and fraction.
         F1 = MOD(D1,1D0)
         F2 = MOD(D2,1D0)
         F = MOD(F1+F2,1D0)
         IF ( F .LT. 0D0 ) F = F+1D0
         D = ANINT(D1-F1) + ANINT(D2-F2) + ANINT(F1+F2-F)
         JD = NINT(D) + 1

!     Express day in Gregorian calendar.
         L = JD + 68569
         N = ( 4*L ) / 146097
         L = L - ( 146097*N + 3 ) / 4
         I = ( 4000 * (L+1) ) / 1461001
         L = L - ( 1461*I ) / 4 + 31
         J = ( 80*L ) / 2447
         ID = L - ( 2447*J ) / 80
         L = J / 11
         IM = J + 2 - 12*L
         IY = 100 * ( N-49 ) + I + L

         FD = F
         J = 0
      END IF

!  Finished.


      END

      DOUBLE PRECISION FUNCTION dau_GMST82 ( DJ1, DJ2 )

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2

      DOUBLE PRECISION DS2R
      PARAMETER ( DS2R = 7.272205216643039903848712D-5 )

!  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

!  Seconds per day, days per Julian century
      DOUBLE PRECISION DAYSEC, CENDAY
      PARAMETER ( DAYSEC = 86400D0, CENDAY = 36525D0 )

!  Coefficients of IAU 1982 GMST-UT1 model
      DOUBLE PRECISION A, B, C, D
      PARAMETER ( A = 24110.54841D0 - DAYSEC/2D0, B = 8640184.812866D0, C = 0.093104D0, D = -6.2D-6 )

!  Note: the first constant, A, has to be adjusted by 12 hours because
!  the UT1 is supplied as a Julian date, which begins at noon.

      DOUBLE PRECISION D1, D2, T, F

      DOUBLE PRECISION iau_ANP

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Julian centuries since fundamental epoch.
      IF ( DJ1 .LT. DJ2 ) THEN
         D1 = DJ1
         D2 = DJ2
      ELSE
         D1 = DJ2
         D2 = DJ1
      END IF
      T = ( D1 + ( D2-DJ0 ) ) / CENDAY

!  Fractional part of JD(UT1), in seconds.
      F = DAYSEC * ( MOD(D1,1D0) + MOD(D2,1D0) )

!  GMST at this UT1.
      dau_GMST82 = iau_ANP ( DS2R * ( (A+(B+(C+D*T)*T)*T) + F ) )

!  Finished.


      END


      DOUBLE PRECISION FUNCTION iau_ANP ( A )

      IMPLICIT NONE

      DOUBLE PRECISION A

!  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      DOUBLE PRECISION W

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      W = MOD(A,D2PI)
      IF ( W .LT. 0D0 ) W = W + D2PI
      iau_ANP = W

!  Finished.


      END

EOF
#=============================================================
#Trabalha com as informações restantes para modificar o arquivo padrão gerado 
    IN=$filename
    str1=${IN#c*_}
    str2=${str1%.dat}
    #echo "Arquivo Original: " $filename " Arquivo Alterado: " $str2
    sed -i 's/arquivo_velho/'$filename'/' PRAIA-results_to_AudeLa.f90
    sed -i 's/OBJETOESTUDADO/'$str2'/' PRAIA-results_to_AudeLa.f90
#=============================================================
#Compila e executa o arquivo padrão gerado
    gfortran PRAIA-results_to_AudeLa.f90 -o PRAIA-results_to_AudeLa
    ./PRAIA-results_to_AudeLa

#=============================================================
# Trabalha com os pontos selecionados em x para modificar o fit
linhai=$(sed '1!d' 'tics_number_gnuplot_corrigido.txt' | xargs) 
linhaii=$(sed '2!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhaiii=$(sed '3!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhaiv=$(sed '4!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhav=$(sed '5!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhavi=$(sed '6!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
#echo "linhas que deverei pegar........" $linhai "," $linhaii "," $linhaiii "," $linhavi "," $linhav "," $linhavi

#=============================================================
# Criar um fortran que seleciona o tempo filtrando quanto tempo passou a mais e depois recebe um tratamento especial para 
# inserir os valores das colunas que serão modificados na sequencia.

cat <<EOF > 'ajustedeescalax.f90'

program comparaarray

	Implicit None
	Real(Kind = 8) :: timemicrosegundos, timemicrosegundosdelta, deltatime
	Integer :: um, dois, tres, quatro, cinco, seis, contalinhasii, j, colunai
	Real(Kind = 8), dimension(100000) :: colunaii, colunaiii, colunatimei, colunatimeii, colunatimeiii
	CHARACTER(20) :: colunaitime, colunaiitime 
	CHARACTER(70) :: colunaiiitime
	CHARACTER(25) :: timeinicial

	OPEN(Unit=27, file="contalinhas.txt", status="Old")
		read(27,*) contalinhasii 
	close (Unit=27)
		um=1
		dois=nint((1.0*(contalinhasii-um)/5.0)+1)
		tres=nint((2.0*(contalinhasii-um)/5.0)+1)
		quatro=nint((3.0*(contalinhasii-um)/5.0)+1)
		cinco=nint((4.0*(contalinhasii-um)/5.0)+1)
		seis=contalinhasii
!		write(*,*) um, dois, tres, quatro, cinco, seis 


	OPEN(Unit=28, file="RESULT_.freemat.00000.dat", status="Old")
	OPEN(Unit=29, file="RESULT_.freemat.00000_2.dat", status="Replace")
	OPEN(Unit=30, file="RESULT_deltatime.dat", status="Old") !deverá ser o arquivo com as diferenças de horários em relação ao primeiro já calculados.

		do j=1,contalinhasii		
			read(28,*) colunai, colunaii(j), colunaiii(j)
			if (colunai == um) then
				read(30,*) colunatimei(j), colunatimeii(j), colunatimeiii(j)
				!write(*,'(I5.1,I5.1,F10.4)') int(colunatimei(j)), int(colunatimeii(j)), real(colunatimeiii(j))
				write(colunaitime,'(I5.1)') int(colunatimei(j))
				write(colunaiitime,'(I5.1)') int(colunatimeii(j))
				write(colunaiiitime,'(F7.4)') real(colunatimeiii(j))
				!write(*,*) colunaitime, colunaiitime, colunaiiitime
				timeinicial = trim(adjustl(colunaitime))//'h'//trim(adjustl(colunaiitime))//'m'//trim(adjustl(colunaiiitime)) //'s'
				timemicrosegundos = real(colunatimei(j))*3600.00 + real(colunatimeii(j))*60.00 + real(colunatimeiii(j))*1.00				
				!write(*,*) timeinicial 
!				write(*,'(A28,F11.4)') "Horário do dia em segundos: ", timemicrosegundos				
				!write(*,'(F11.4)') timemicrosegundos
				write(29,10043) colunai, colunaii(j), colunaiii(j), timeinicial
			else if (colunai == dois .OR. colunai == tres & 
.OR. colunai == quatro .OR. colunai == cinco .OR. colunai == seis) then
				read(30,*) colunatimei(j), colunatimeii(j), colunatimeiii(j)
				!write(*,'(I5.1,I5.1,F10.4)') int(colunatimei(j)), int(colunatimeii(j)), real(colunatimeiii(j))
				write(colunaitime,'(I5.1)') int(colunatimei(j))
				write(colunaiitime,'(I5.1)') int(colunatimeii(j))
				write(colunaiiitime,'(F7.4)') real(colunatimeiii(j))
				!write(*,*) colunaitime, colunaiitime, colunaiiitime
				timeinicial = trim(adjustl(colunaitime))//'h'//trim(adjustl(colunaiitime))//'m'//trim(adjustl(colunaiiitime)) //'s'
				timemicrosegundosdelta = real(colunatimei(j))*3600.00 + real(colunatimeii(j))*60.00 + real(colunatimeiii(j))*1.00			
				deltatime = timemicrosegundosdelta - timemicrosegundos				
!				write(*,*) timeinicial
				!write(*,'(A24,F11.4)') "Incremento (segundos): ", deltatime
				write(29,10044) colunai, colunaii(j), colunaiii(j), deltatime
			else	
				read(30,*) colunatimei(j), colunatimeii(j), colunatimeiii(j)
				write(29,10045) colunai, colunaii(j), colunaiii(j)
			end if

		end do
	close (Unit=28)
	close (Unit=29)
	close (Unit=30)

10043 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4 , ' ', A25)
10044 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4 , ' ', F15.4 ,'s')
10045 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4)


end program comparaarray

EOF

	gfortran ajustedeescalax.f90 -o ajustedeescalax
	./ajustedeescalax

#=============================================================
#Chama o Gnuplot e manda plotar os Gráficos

cat <<EOF >'plot_'$str2'.gnu'
set grid #Insere Grades
set key box #Legenda dentro de um quadro
#set key out vert #Legenda Fora do Grafico
set key bottom right #Legenda na parte inferior a direita
set title "Prévia de Resultados da Análise da\nCurva de Luz do Corpo Celeste: $str2"
set xlabel "Horário inicial (UTC) com Incrementos em segundos"
set ylabel "Fluxo de Luz do Objeto"
set terminal pngcairo size 700,524 enhanced font 'Verdana,8'
set lmargin 12
set rmargin 8


set output "RESULT_plot_$str2.png"

#set xtics rotate by 330
#set xtics ("$linhai" $linhai,"+$linhaii" $linhaii,"+$linhaiii" $linhaiii,"+$linhaiv" $linhaiv,"+$linhav" $linhav,"+$linhavi" $linhavi)

plot 'RESULT_.freemat.00000_2.dat' using 1:2:xtic(4) title 'Curva de Luz' with linespoints ls 3 pt 2
#, 'RESULT_.freemat.00000.dat' using 1:3 title 'Star' with linespoints ls 4
#, 'RESULT_.freemat.00000.dat' using 1:2 smooth bezier notitle ls 1

EOF

   # manda gnuplot 
   gnuplot 'plot_'$str2'.gnu'
#=============================================================
#Remove o arquivo .f90, o gerado executavel e demais arquivos inuteis para aquela pasta.
   rm -f PRAIA-results_to_AudeLa.f90 PRAIA-results_to_AudeLa 'plot_'$str2'.gnu' contalinhas.txt intervalo_de_tempo_gnuplot.txt tics_number_gnuplot.txt arrendadomento_de_quadros.f90 arrendadomento_de_quadros tics_number_gnuplot.txt tics_number_gnuplot_corrigido.txt ajustedeescalax.f90 ajustedeescalax RESULT_.freemat.00000_2.dat RESULT_deltatime.dat

   # abre o arquivo do gráfico plotado
   shotwell "RESULT_plot_$str2.png" &
#=============================================================
done









#=====================================================================
#Se foi escolhido a opção de fotometria (2)
else
	if (("$opcaodetrabalho" == "2"))
	then
#=====================================================================
#Procura quantas linhas dará
gawk '{ sum += $1 }; END { print sum }' photometry_*.dat > contalinhas.txt

#Cria Fit File
gawk '{ sum2 += $1 }; END { print (sum2 / 5)}' photometry_*.dat > intervalo_de_tempo_gnuplot.txt
contouaslinhas=$(<intervalo_de_tempo_gnuplot.txt)
echo "1" > tics_number_gnuplot.txt
echo "$contouaslinhas * 1" | bc -l >> tics_number_gnuplot.txt
echo "$contouaslinhas * 2" | bc -l >> tics_number_gnuplot.txt
echo "$contouaslinhas * 3" | bc -l >> tics_number_gnuplot.txt
echo "$contouaslinhas * 4" | bc -l >> tics_number_gnuplot.txt
gawk '{ sum += $1 }; END { print sum }' photometry_*.dat >> tics_number_gnuplot.txt

#Arredonda os valores
cat <<EOF >arrendadomento_de_quadros.f90
Program Arredondamento
	Implicit None
! ======================================================================================
	Integer :: contaarredondador	
	Real(Kind = 8), dimension(6) :: arredondar
	OPEN(Unit=21, file="tics_number_gnuplot.txt", status="Old")
	OPEN(Unit=22, file="tics_number_gnuplot_corrigido.txt", status="Replace")
	do contaarredondador=1, 6
		read(21,*) arredondar(contaarredondador)
		Write(22,*) ceiling(arredondar(contaarredondador))		 
	end do
	close (Unit=22)	
	close (Unit=21)
End Program Arredondamento
EOF

#Encerra Fit File
	gfortran 'arrendadomento_de_quadros.f90' -o 'arrendadomento_de_quadros'
	./arrendadomento_de_quadros


for filename in $(ls -c photometry_*.dat | head -1);
do
#=============================================================
#Insere o seguinte código-fonte em um arquivo do Fortran Gerado como Padrão
cat <<EOF >PRAIA-results_to_AudeLa.f90
Program PRAIAtoAudeLA
	Implicit None
! ======================================================================================
	Real(Kind = 8) :: U02, U03, U04, U05, U07, U08, U10, U11, U13, mediacalibradores
	DOUBLE PRECISION, dimension(100000) :: JDay
	Real(Kind = 8), dimension(100000) :: BrilhoObjeto , BrilhoRefI, BrilhoRefII, BrilhoRefIII 
	Integer :: i, j, a, b, contalinhas, mediaverificada
	CHARACTER(20) :: NomeDoFit
	CHARACTER(13) :: TimeISOa
	CHARACTER(2) :: TimeISOb
	CHARACTER(7) :: TimeISOc
	CHARACTER(27) :: TimeISO
	 
! ======================================================================================
!Arquivo: arquivo_velho
TimeISO = ''
mediaverificada=0
! Leitura nos arquivos de dados
	OPEN(Unit=15, file="contalinhas.txt", status="Old")
		read(15,*) contalinhas 
	close (Unit=15)

	OPEN(Unit=2, file="photometry_OBJETOESTUDADO.dat", status="Old")
	! Gerando arquivo corrigido
	  OPEN(Unit=3, file="RESULT_OBJETOESTUDADO_.00000.csv", status="Replace")
	  OPEN(Unit=99, file="RESULT_OBJETOESTUDADO_.time", status="Replace")
	  OPEN(Unit=98, file="RESULT_.freemat.00000.dat", status="Replace")
	  OPEN(Unit=97, file="RESULT_deltatime.dat", status="Replace")
		Write(3,*) "# ** atos - Audela - Linux  * "
		Write(99,*) "# ** atos - Audela - Linux  * "
		Write(3,*) "# FPS 25"
		Write(99,*) "# FPS 25"		
		Write(3,*) "idframe,jd,dateiso,obj_fint     ,obj_pixmax   ,obj_intensite, &
obj_sigmafond,obj_snint    ,obj_snpx     ,obj_delta    ,obj_xpos,obj_ypos,obj_xfwhm,obj_yfwhm,&
ref_fint     ,ref_pixmax   ,ref_intensite,ref_sigmafond,ref_snint    ,ref_snpx     ,ref_delta    ,&
ref_xpos,ref_ypos,ref_xfwhm,ref_yfwhm,img_intmin ,img_intmax,img_intmoy,img_sigma ,img_xsize,img_ysize"
		Write(99,*) "idframe, jd, dateiso, verif, ocr, interpol"
		do i=1, contalinhas
	  		read(2,*) a, b, JDay(i), U02, BrilhoObjeto(i), U04, U05, BrilhoRefI(i), U07, U08, BrilhoRefII(i), U10,U11, BrilhoRefIII(i), U13, U03, NomeDoFit
			if(mediaverificada < 1) then
				write (*,*) "Estrelas/Objetos Calibradores Encontrados", b
				write (*,*) "Deseja utilizar Utilizar Quantos Objetos Calibradores? (Exemplo: 2)"				
				read (*,*) mediacalibradores
				mediaverificada = 2.0
			end if
			if(mediacalibradores == 1.0 .OR. mediacalibradores < 1.0 .OR. mediacalibradores > 3.0) then
				BrilhoRefI(i) = BrilhoRefI(i)
			else if(mediacalibradores == 2.0) then
				BrilhoRefI(i) = ((BrilhoRefI(i) + BrilhoRefII(i))/2.0)
			else if(mediacalibradores == 3.0) then
				BrilhoRefI(i) = ((BrilhoRefI(i) + BrilhoRefII(i) + BrilhoRefIII(i))/3.0)
			end if

			Call jd2greg(JDay(i))
			
			OPEN(Unit=11, file="auxiliar2.dat", status="Old")
				do j=1,1
					read(11,*) TimeISOa, TimeISOb, TimeISOc
					!concatenação
					TimeISO = TimeISOa //':'// TimeISOb //':'// TimeISOc
				end do
			close (Unit=11)

			write(97,*) TimeISO(12:13), "    " , TimeISO(15:16), "    " , TimeISO(18:24)

			write(3,10010) i, JDay(i), TimeISO, BrilhoObjeto(i), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, BrilhoRefI(i), 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
			write(99,10030) i, JDay(i), TimeISO, 1, 1, 1

			write(98,10020) i, BrilhoObjeto(i), BrilhoRefI(i)
	   	end do
	close(97)	
	close(98)	
	close(99)	  
	close(3)
	close(Unit = 2)



	!Depois de tudo feito
	OPEN(Unit=11, file="auxiliar.dat", status="Old")
	close (Unit=11, status='delete')
	OPEN(Unit=11, file="auxiliar2.dat", status="Old")
	close (Unit=11, status='delete')

	Write(*,*) " "
	Write(*,*) "Resultados de Fotometria foram salvos com sucesso!"
	Write(*,*) "=========================================================================================="

	Write(*,*) " "
	Write(*,*) " "

! ======================================================================================
! Formatos de Saida
10010 FORMAT(I5.1, ',', F16.8, ',', A25, ',', F15.4, ',', I3.1, &
',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', F15.4, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', I3.1, ',', &
I3.1, ',', I3.1, ',', I3.1)

10020 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4)

10030 FORMAT(I5.1, ',', F16.8, ',', A25, ',', I3.1, ',', I3.1, ',', I3.1)

10000 FORMAT(' ', I5.1, ',', F16.8, ',', A25, ',', F10.4, ',', I3.1, ',', A30) !Relatorio Final
10001 FORMAT('i = ',I3.1) !para inteiro
10002 FORMAT(5x,'JD(i) = ',F16.8) !para float de JD
10003 FORMAT(5x,'BrilhoObjeto(i) = ',F9.4) !float de brilho
10004 FORMAT(5x,'Nome do Fit ='A30) !texto

! ======================================================================================

End Program PRAIAtoAudeLA
! ======================================================================================
! Subrotinas que convertem uma data JD para ISO (Julian Day to a Gregorian Date)

SUBROUTINE jd2greg (julianTime)
	!DOUBLE PRECISION :: julianTime
	!julianTime = 2*julianTime
	!write(*,20002) julianTime
	!20002 FORMAT(5x,'JDay(i) = ',F16.8) !para float de JD

      IMPLICIT REAL *8 (A-H,O-Z)
      DOUBLE PRECISION :: julianTime
      character*24 ler, julianTimeText
      hmsgms(ai,aj,aa)=ai+aj/60.d0+aa/3600.d0
      dj2000=2451544.5d0
      dj1=2400000.5D0
      idim=250

      julianTimeText = ''
    open(4,file="auxiliar.dat")
       write(4,*) julianTime
    close(4)

    open(4,file="auxiliar.dat")
       read(4, *) julianTimeText
       !write(*,*) A_char
    close(4)

      ler=julianTimeText !2456977.23464133'
!'
      read (ler,*) data
      iflag=1
      if (data.lt.2500.d0) iflag=2
      if (data.gt.25.d0.and.data.lt.2500.d0) iflag=3
      if (iflag.ne.1) go to 30
 15   continue
      djm=data-dj1
      call iau_jd2cal (dj1,djm,iutano,iutmes,iutdia,fd,jjj)
      hora=fd*24.d0
      iuth=hora
      iutm=(hora-iuth)*60.d0
      sut=((hora-iuth)*60.d0-iutm)*60.d0
      frano=(data-dj2000)/365.25d0+2000.d0
      ler=''
!      write (ler,20) iuth,iutm,sut,iutdia,iutmes,iutano,data,frano
! 20   format(2(i2.2,':'),f7.4,':',2(i2.2,':'),i4.4,':',f18.10,':',
!     ?f15.10)

!      write (ler,20) iuth,iutm,sut,iutdia,iutmes,iutano
    open(666,file="auxiliar2.dat")
 !      write (666,30001) iutano,iutmes,iutdia,iuth,iutm,sut
 !30001   format(i4.4,'-',1(i2.2,'-'),i2.2,'T', 2(i2.2,':'),f7.4)
   ! close(666)


      write (ler,20) iutano,iutmes,iutdia,iuth,iutm,sut
!	write (*,20) iutano,iutmes,iutdia,iuth,iutm,sut
 20   format(i4.4,'-',1(i2.2,'-'),i2.2,'T', 2(i2.2,':'),f7.4)
!20   format(i4.4,'-',1(i2.2,'-'),i2.2,'T', 2(i2.2,':'),f7.4)


      do i=idim,1,-1
      if (ler(i:i).ne.' ') go to 25
      enddo
 25   ii=i

      do i=1,ii-1
      if (ler(i:i).eq.' ') ler(i:i)='0'
      enddo

      do i=1,ii-1
      if (ler(i:i).eq.':') ler(i:i)=' '
      enddo

!      write (*,*) ler
	write (666,*) ler
	close(666)
      go to 100

 30   continue

      if (iflag.eq.3) then
      data=365.25d0*(data-2000.d0)+dj2000
      go to 15
      endif

      read (ler,*) uth,utm,sut,utdia,iutmes,iutano

      fd=hmsgms(uth,utm,sut)/24.d0

      iutdia=utdia

      fdd=utdia-iutdia

      fd=fd+fdd

      call iau_CAL2JD (iutano,iutmes,iutdia,djm0,djm,jjj)


      djm=djm+fd
            
      data=djm+djm0


      go to 15

 100  continue




END SUBROUTINE jd2greg



      SUBROUTINE iau_CAL2JD ( IY, IM, ID, DJM0, DJM, J )

      IMPLICIT NONE

      INTEGER IY, IM, ID
      DOUBLE PRECISION DJM0, DJM
      INTEGER J, MY, IYPMY

!  Earliest year allowed (4800BC)
      INTEGER IYMIN
      PARAMETER ( IYMIN = -4799 )

!  Month lengths in days
      INTEGER MTAB(12)
      DATA MTAB / 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Preset status.
      J = 0

!  Validate year.
      IF ( IY.LT.IYMIN ) THEN
         J = -1
      ELSE

!     Validate month.
         IF ( IM.GE.1 .AND. IM.LE.12 ) THEN

!        Allow for leap year.
            IF ( MOD(IY,4) .EQ. 0 ) THEN
               MTAB(2) = 29
            ELSE
               MTAB(2) = 28
            END IF
            IF ( MOD(IY,100).EQ.0 .AND. MOD(IY,400).NE.0 ) MTAB(2) = 28

!        Validate day.
            IF ( ID.LT.1 .OR. ID.GT.MTAB(IM) ) J = -3

!        Result.
            MY = ( IM - 14 ) / 12
            IYPMY = IY + MY
            DJM0 = 2400000.5D0
            DJM = DBLE( (1461*(IYPMY+4800))/4 + (367*(IM-2-(12*MY)))/12 - (3*((IYPMY+4900)/100))/4 + ID-2432076)

!        Bad month
         ELSE
            J = -2
         END IF
      END IF


      END


      SUBROUTINE iau_jd2cal ( DJ1, DJ2, IY, IM, ID, FD, J )

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2
      INTEGER IY, IM, ID
      DOUBLE PRECISION FD
      INTEGER J

!  Minimum and maximum allowed JD
      DOUBLE PRECISION DJMIN, DJMAX
      PARAMETER ( DJMIN = -68569.5D0, DJMAX = 1D9 )

      INTEGER JD, L, N, I
      DOUBLE PRECISION DJ, D1, D2, F1, F2, F, D

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Check if date is acceptable.
      DJ = DJ1 + DJ2
      IF ( DJ.LT.DJMIN .OR. DJ.GT.DJMAX ) THEN
         J = -1
      ELSE
         J = 0

!     Copy the date, big then small, and re-align to midnight.
         IF ( DJ1 .GE. DJ2 ) THEN
            D1 = DJ1
            D2 = DJ2
         ELSE
            D1 = DJ2
            D2 = DJ1
         END IF
         D2 = D2 - 0.5D0

!     Separate day and fraction.
         F1 = MOD(D1,1D0)
         F2 = MOD(D2,1D0)
         F = MOD(F1+F2,1D0)
         IF ( F .LT. 0D0 ) F = F+1D0
         D = ANINT(D1-F1) + ANINT(D2-F2) + ANINT(F1+F2-F)
         JD = NINT(D) + 1

!     Express day in Gregorian calendar.
         L = JD + 68569
         N = ( 4*L ) / 146097
         L = L - ( 146097*N + 3 ) / 4
         I = ( 4000 * (L+1) ) / 1461001
         L = L - ( 1461*I ) / 4 + 31
         J = ( 80*L ) / 2447
         ID = L - ( 2447*J ) / 80
         L = J / 11
         IM = J + 2 - 12*L
         IY = 100 * ( N-49 ) + I + L

         FD = F
         J = 0
      END IF

!  Finished.


      END

      DOUBLE PRECISION FUNCTION dau_GMST82 ( DJ1, DJ2 )

      IMPLICIT NONE

      DOUBLE PRECISION DJ1, DJ2

      DOUBLE PRECISION DS2R
      PARAMETER ( DS2R = 7.272205216643039903848712D-5 )

!  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

!  Seconds per day, days per Julian century
      DOUBLE PRECISION DAYSEC, CENDAY
      PARAMETER ( DAYSEC = 86400D0, CENDAY = 36525D0 )

!  Coefficients of IAU 1982 GMST-UT1 model
      DOUBLE PRECISION A, B, C, D
      PARAMETER ( A = 24110.54841D0 - DAYSEC/2D0, B = 8640184.812866D0, C = 0.093104D0, D = -6.2D-6 )

!  Note: the first constant, A, has to be adjusted by 12 hours because
!  the UT1 is supplied as a Julian date, which begins at noon.

      DOUBLE PRECISION D1, D2, T, F

      DOUBLE PRECISION iau_ANP

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Julian centuries since fundamental epoch.
      IF ( DJ1 .LT. DJ2 ) THEN
         D1 = DJ1
         D2 = DJ2
      ELSE
         D1 = DJ2
         D2 = DJ1
      END IF
      T = ( D1 + ( D2-DJ0 ) ) / CENDAY

!  Fractional part of JD(UT1), in seconds.
      F = DAYSEC * ( MOD(D1,1D0) + MOD(D2,1D0) )

!  GMST at this UT1.
      dau_GMST82 = iau_ANP ( DS2R * ( (A+(B+(C+D*T)*T)*T) + F ) )

!  Finished.


      END


      DOUBLE PRECISION FUNCTION iau_ANP ( A )

      IMPLICIT NONE

      DOUBLE PRECISION A

!  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

      DOUBLE PRECISION W

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      W = MOD(A,D2PI)
      IF ( W .LT. 0D0 ) W = W + D2PI
      iau_ANP = W

!  Finished.


      END

EOF
#=============================================================
#Trabalha com as informações restantes para modificar o arquivo padrão gerado 
    IN=$filename
    str1=${IN#p*_}
    str2=${str1%.dat}
#    echo "Arquivo Original: " $filename " Arquivo Alterado: " $str2
    sed -i 's/arquivo_velho/'$filename'/' PRAIA-results_to_AudeLa.f90
    sed -i 's/OBJETOESTUDADO/'$str2'/' PRAIA-results_to_AudeLa.f90
#=============================================================
#Compila e executa o arquivo padrão gerado
    gfortran PRAIA-results_to_AudeLa.f90 -o PRAIA-results_to_AudeLa
    ./PRAIA-results_to_AudeLa

#=============================================================
# Trabalha com os pontos selecionados em x para modificar o fit
linhai=$(sed '1!d' 'tics_number_gnuplot_corrigido.txt' | xargs) 
linhaii=$(sed '2!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhaiii=$(sed '3!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhaiv=$(sed '4!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhav=$(sed '5!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
linhavi=$(sed '6!d' 'tics_number_gnuplot_corrigido.txt'| xargs) 
#echo "linhas que deverei pegar........" $linhai "," $linhaii "," $linhaiii "," $linhavi "," $linhav "," $linhavi

#=============================================================
# Criar um fortran que seleciona o tempo filtrando quanto tempo passou a mais e depois recebe um tratamento especial para 
# inserir os valores das colunas que serão modificados na sequencia.

cat <<EOF > 'ajustedeescalax.f90'

program comparaarray

	Implicit None
	Real(Kind = 8) :: timemicrosegundos, timemicrosegundosdelta, deltatime
	Integer :: um, dois, tres, quatro, cinco, seis, contalinhasii, j, colunai
	Real(Kind = 8), dimension(100000) :: colunaii, colunaiii, colunatimei, colunatimeii, colunatimeiii
	CHARACTER(20) :: colunaitime, colunaiitime 
	CHARACTER(70) :: colunaiiitime
	CHARACTER(25) :: timeinicial

	OPEN(Unit=27, file="contalinhas.txt", status="Old")
		read(27,*) contalinhasii 
	close (Unit=27)
		um=1
		dois=nint((1.0*(contalinhasii-um)/5.0)+1)
		tres=nint((2.0*(contalinhasii-um)/5.0)+1)
		quatro=nint((3.0*(contalinhasii-um)/5.0)+1)
		cinco=nint((4.0*(contalinhasii-um)/5.0)+1)
		seis=contalinhasii
!		write(*,*) um, dois, tres, quatro, cinco, seis 


	OPEN(Unit=28, file="RESULT_.freemat.00000.dat", status="Old")
	OPEN(Unit=29, file="RESULT_.freemat.00000_2.dat", status="Replace")
	OPEN(Unit=30, file="RESULT_deltatime.dat", status="Old") !deverá ser o arquivo com as diferenças de horários em relação ao primeiro já calculados.

		do j=1,contalinhasii		
			read(28,*) colunai, colunaii(j), colunaiii(j)
			if (colunai == um) then
				read(30,*) colunatimei(j), colunatimeii(j), colunatimeiii(j)
				!write(*,'(I5.1,I5.1,F10.4)') int(colunatimei(j)), int(colunatimeii(j)), real(colunatimeiii(j))
				write(colunaitime,'(I5.1)') int(colunatimei(j))
				write(colunaiitime,'(I5.1)') int(colunatimeii(j))
				write(colunaiiitime,'(F7.4)') real(colunatimeiii(j))
				!write(*,*) colunaitime, colunaiitime, colunaiiitime
				timeinicial = trim(adjustl(colunaitime))//'h'//trim(adjustl(colunaiitime))//'m'//trim(adjustl(colunaiiitime)) //'s'
				timemicrosegundos = real(colunatimei(j))*3600.00 + real(colunatimeii(j))*60.00 + real(colunatimeiii(j))*1.00				
				!write(*,*) timeinicial 
!				write(*,'(A28,F11.4)') "Horário do dia em segundos: ", timemicrosegundos				
				!write(*,'(F11.4)') timemicrosegundos
				write(29,10043) colunai, colunaii(j), colunaiii(j), timeinicial
			else if (colunai == dois .OR. colunai == tres & 
.OR. colunai == quatro .OR. colunai == cinco .OR. colunai == seis) then
				read(30,*) colunatimei(j), colunatimeii(j), colunatimeiii(j)
				!write(*,'(I5.1,I5.1,F10.4)') int(colunatimei(j)), int(colunatimeii(j)), real(colunatimeiii(j))
				write(colunaitime,'(I5.1)') int(colunatimei(j))
				write(colunaiitime,'(I5.1)') int(colunatimeii(j))
				write(colunaiiitime,'(F7.4)') real(colunatimeiii(j))
				!write(*,*) colunaitime, colunaiitime, colunaiiitime
				timeinicial = trim(adjustl(colunaitime))//'h'//trim(adjustl(colunaiitime))//'m'//trim(adjustl(colunaiiitime)) //'s'
				timemicrosegundosdelta = real(colunatimei(j))*3600.00 + real(colunatimeii(j))*60.00 + real(colunatimeiii(j))*1.00			
				deltatime = timemicrosegundosdelta - timemicrosegundos				
!				write(*,*) timeinicial
				!write(*,'(A24,F11.4)') "Incremento (segundos): ", deltatime
				write(29,10044) colunai, colunaii(j), colunaiii(j), deltatime
			else	
				read(30,*) colunatimei(j), colunatimeii(j), colunatimeiii(j)
				write(29,10045) colunai, colunaii(j), colunaiii(j)
			end if

		end do
	close (Unit=28)
	close (Unit=29)
	close (Unit=30)

10043 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4 , ' ', A25)
10044 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4 , ' ', F15.4 ,'s')
10045 FORMAT(' ', I5.1, ' ', F15.4, ' ', F15.4)


end program comparaarray

EOF

	gfortran ajustedeescalax.f90 -o ajustedeescalax
	./ajustedeescalax

#=============================================================
#Chama o Gnuplot e manda plotar os Gráficos

cat <<EOF >'plot_'$str2'.gnu'
set grid #Insere Grades
set key box #Legenda dentro de um quadro
#set key out vert #Legenda Fora do Grafico
set key bottom right #Legenda na parte inferior a direita
set title "Prévia de Resultados da Análise da\nCurva de Luz do Corpo Celeste: $str2"
set xlabel "Horário inicial (UTC) com Incrementos em segundos"
set ylabel "Fluxo de Luz do Objeto"
set terminal pngcairo size 700,524 enhanced font 'Verdana,8'
set lmargin 12
set rmargin 8


set output "RESULT_plot_$str2.png"

#set xtics rotate by 330
#set xtics ("$linhai" $linhai,"+$linhaii" $linhaii,"+$linhaiii" $linhaiii,"+$linhaiv" $linhaiv,"+$linhav" $linhav,"+$linhavi" $linhavi)

plot 'RESULT_.freemat.00000_2.dat' using 1:2:xtic(4) title 'Curva de Luz' with linespoints ls 3 pt 2
#, 'RESULT_.freemat.00000.dat' using 1:3 title 'Star' with linespoints ls 4
#, 'RESULT_.freemat.00000.dat' using 1:2 smooth bezier notitle ls 1

EOF

   # manda gnuplot 
   gnuplot 'plot_'$str2'.gnu'
#=============================================================
#Remove o arquivo .f90, o gerado executavel e demais arquivos inuteis para aquela pasta.
   rm -f PRAIA-results_to_AudeLa.f90 PRAIA-results_to_AudeLa 'plot_'$str2'.gnu' contalinhas.txt intervalo_de_tempo_gnuplot.txt tics_number_gnuplot.txt arrendadomento_de_quadros.f90 arrendadomento_de_quadros tics_number_gnuplot.txt tics_number_gnuplot_corrigido.txt ajustedeescalax.f90 ajustedeescalax RESULT_.freemat.00000_2.dat RESULT_deltatime.dat

   # abre o arquivo do gráfico plotado
   shotwell "RESULT_plot_$str2.png" &
#=============================================================
done

#=====================================================================
	fi
fi
