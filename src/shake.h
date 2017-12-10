c global variables for Shake routine
c
        integer iShake           ! shake option
	integer nbondShake
        integer shitMX           ! MAX numb of shake iteration
        real shTol               ! accuracy
        integer shExitFlag       ! exit flag
c
        common/shake01/iShake,nbondShake,shitMX,shTol,shExitFlag
c
