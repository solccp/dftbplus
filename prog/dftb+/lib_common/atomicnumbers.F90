module atomicnumbers
  use accuracy
  
  public
  type atomicnumberpair
     character(2) :: symbol
     integer      :: atomicnumber
  end type atomicnumberpair

  integer, parameter :: nAtomicNumberPair = 106;


  type(atomicnumberpair), parameter :: atomicNumberPairs(nAtomicNumberPair) = (/ &
    &atomicnumberpair("X ", 0), & 
    &atomicnumberpair("H ", 1), & 
    &atomicnumberpair("D ", 1), & 
    &atomicnumberpair("T ", 1), & 
    &atomicnumberpair("He", 2), & 
    &atomicnumberpair("Li", 3), & 
    &atomicnumberpair("Be", 4), & 
    &atomicnumberpair("B ", 5), & 
    &atomicnumberpair("C ", 6), & 
    &atomicnumberpair("N ", 7), & 
    &atomicnumberpair("O ", 8), & 
    &atomicnumberpair("F ", 9), & 
    &atomicnumberpair("Ne", 10), & 
    &atomicnumberpair("Na", 11), & 
    &atomicnumberpair("Mg", 12), & 
    &atomicnumberpair("Al", 13), & 
    &atomicnumberpair("Si", 14), & 
    &atomicnumberpair("P ", 15), & 
    &atomicnumberpair("S ", 16), & 
    &atomicnumberpair("Cl", 17), & 
    &atomicnumberpair("Ar", 18), & 
    &atomicnumberpair("K ", 19), & 
    &atomicnumberpair("Ca", 20), & 
    &atomicnumberpair("Sc", 21), & 
    &atomicnumberpair("Ti", 22), & 
    &atomicnumberpair("V ", 23), & 
    &atomicnumberpair("Cr", 24), & 
    &atomicnumberpair("Mn", 25), & 
    &atomicnumberpair("Fe", 26), & 
    &atomicnumberpair("Co", 27), & 
    &atomicnumberpair("Ni", 28), & 
    &atomicnumberpair("Cu", 29), & 
    &atomicnumberpair("Zn", 30), & 
    &atomicnumberpair("Ga", 31), & 
    &atomicnumberpair("Ge", 32), & 
    &atomicnumberpair("As", 33), & 
    &atomicnumberpair("Se", 34), & 
    &atomicnumberpair("Br", 35), & 
    &atomicnumberpair("Kr", 36), & 
    &atomicnumberpair("Rb", 37), & 
    &atomicnumberpair("Sr", 38), & 
    &atomicnumberpair("Y ", 39), & 
    &atomicnumberpair("Zr", 40), & 
    &atomicnumberpair("Nb", 41), & 
    &atomicnumberpair("Mo", 42), & 
    &atomicnumberpair("Tc", 43), & 
    &atomicnumberpair("Ru", 44), & 
    &atomicnumberpair("Rh", 45), & 
    &atomicnumberpair("Pd", 46), & 
    &atomicnumberpair("Ag", 47), & 
    &atomicnumberpair("Cd", 48), & 
    &atomicnumberpair("In", 49), & 
    &atomicnumberpair("Sn", 50), & 
    &atomicnumberpair("Sb", 51), & 
    &atomicnumberpair("Te", 52), & 
    &atomicnumberpair("I ", 53), & 
    &atomicnumberpair("Xe", 54), & 
    &atomicnumberpair("Cs", 55), & 
    &atomicnumberpair("Ba", 56), & 
    &atomicnumberpair("La", 57), & 
    &atomicnumberpair("Ce", 58), & 
    &atomicnumberpair("Pr", 59), & 
    &atomicnumberpair("Nd", 60), & 
    &atomicnumberpair("Pm", 61), & 
    &atomicnumberpair("Sm", 62), & 
    &atomicnumberpair("Eu", 63), & 
    &atomicnumberpair("Gd", 64), & 
    &atomicnumberpair("Tb", 65), & 
    &atomicnumberpair("Dy", 66), & 
    &atomicnumberpair("Ho", 67), & 
    &atomicnumberpair("Er", 68), & 
    &atomicnumberpair("Tm", 69), & 
    &atomicnumberpair("Yb", 70), & 
    &atomicnumberpair("Lu", 71), & 
    &atomicnumberpair("Hf", 72), & 
    &atomicnumberpair("Ta", 73), & 
    &atomicnumberpair("W ", 74), & 
    &atomicnumberpair("Re", 75), & 
    &atomicnumberpair("Os", 76), & 
    &atomicnumberpair("Ir", 77), & 
    &atomicnumberpair("Pt", 78), & 
    &atomicnumberpair("Au", 79), & 
    &atomicnumberpair("Hg", 80), & 
    &atomicnumberpair("Tl", 81), & 
    &atomicnumberpair("Pb", 82), & 
    &atomicnumberpair("Bi", 83), & 
    &atomicnumberpair("Po", 84), & 
    &atomicnumberpair("At", 85), & 
    &atomicnumberpair("Rn", 86), & 
    &atomicnumberpair("Fr", 87), & 
    &atomicnumberpair("Ra", 88), & 
    &atomicnumberpair("Ac", 89), & 
    &atomicnumberpair("Th", 90), & 
    &atomicnumberpair("Pa", 91), & 
    &atomicnumberpair("U ", 92), & 
    &atomicnumberpair("Np", 93), & 
    &atomicnumberpair("Pu", 94), & 
    &atomicnumberpair("Am", 95), & 
    &atomicnumberpair("Cm", 96), & 
    &atomicnumberpair("Bk", 97), & 
    &atomicnumberpair("Cf", 98), & 
    &atomicnumberpair("Es", 99), & 
    &atomicnumberpair("Fm", 100), & 
    &atomicnumberpair("Md", 101), & 
    &atomicnumberpair("No", 102), & 
    &atomicnumberpair("Lr", 103) & 
      &/)

  contains
  function getAtomicNumberFromSymbol(symbol_)
    integer getAtomicNumberFromSymbol
    character(mc), intent(in) :: symbol_
    integer :: ii

    do ii = 1, nAtomicNumberPair
        if (trim(atomicNumberPairs(ii)%symbol) == trim(symbol_) ) then
            getAtomicNumberFromSymbol = atomicNumberPairs(ii)%atomicnumber
            return
        end if
    end do
    getAtomicNumberFromSymbol = 0 
  end function
end module atomicnumbers

