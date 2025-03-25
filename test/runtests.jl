using Test
using Mosek
using MosekTools

########################
# test mosek installation
@test begin
    x = maketask() do task
        # Use remote server: putoptserverhost(task,"http://solve.mosek.com:30080")
        appendvars(task, 1)                             # 1 variable x
        putcj(task, 1, 1.0)                             # c_0 = 1.0
        putvarbound(task, 1, MSK_BK_RA, 2.0, 3.0)       # 2.0 <= x <= 3.0
        putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE) # minimize
    
        optimize(task)                                  # Optimize
    
        getxx(task, MSK_SOL_ITR)                    # Get solution
    end
    x == [2.0]
end