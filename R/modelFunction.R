modelFunction <- function(dose, parm2mat, drcFct, cm, assayNoOld, upperPos, retFct, doseScaling, respScaling, isFinite)
{
    if (!is.null(retFct))
    {
        drcFct <- retFct(doseScaling, respScaling)
    }
    drcFct1 <- function(dose, parm)
    {
        drcFct(dose, (parm2mat(parm))[isFinite, , drop = FALSE])
    }

    if (is.null(cm))
    {
        multCurves <- function(dose, parm)
        {
           drcFct1(dose, parm)
        }
    } else {  # not adapting to scaling (not using drcFct1)!!!
        iv <- isFinite & (assayNoOld == cm)
        niv <- !iv
        fctEval <- rep(0, length(dose))

        multCurves <- function(dose, parm)
        {
            parmVal <- (parm2mat(parm))[isFinite, , drop = FALSE]
            fctEval[iv] <- parmVal[iv, upperPos, drop = FALSE]
            fctEval[niv] <- drcFct(dose[niv], parmVal[niv, , drop = FALSE])

            fctEval
        }
    }
    
    multCurves
}

