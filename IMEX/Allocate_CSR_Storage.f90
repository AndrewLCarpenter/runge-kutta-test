      subroutine Allocate_CSR_Storage(problem,n)

        use precision_vars
        use CSR_Variables

        implicit none

        integer,                        intent(in   )  :: problem
        integer,                        intent(  out)  :: n

        if(problem <= 4) then
            nJac = 2 
          nnzJac = 4

        elseif(problem == 5) then
            nJac = 3 
          nnzJac = 9

        endif
        n = nJac

        allocate(iaJac(nJac+1))
        allocate(jaJac(nnzJac))
        allocate( aJac(nnzJac))

        allocate( juJac(nJac+1))
        allocate(jLUJac(nnzJac))
        allocate(aLUJac(nnzJac))

        allocate(iw(nJac))

      end subroutine
