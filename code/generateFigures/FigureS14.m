function FigureS14()
% Save figures comparing heritabitity in bins

[RvsF_BOTH, FvsP_BOTH] = compare_numberExcludedSubjects('BOTH');
figureName = 'makeFigures/Heritability_subjectBINS.png';
print(gcf,figureName,'-dpng','-r600');


[RvsF_VAR, FvsP_VAR] = compare_numberExcludedSubjects('VARrem');
figureName = 'makeFigures/Heritability_varianceBINS.png';
print(gcf,figureName,'-dpng','-r600');

end