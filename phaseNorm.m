function ComplexFieldNorm = phaseNorm(ComplexField)
    arguments
        ComplexField
    end

    phase = angle(ComplexField);

    MIN = min(min(phase)); % idk why double min

    ComplexFieldNorm = ComplexField.*exp(-1j*MIN);
    
end