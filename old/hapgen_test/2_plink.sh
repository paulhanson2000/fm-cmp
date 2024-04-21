# Combine the case/control HAPGEN output files.
plink --data out.cases    --oxford-single-chr 1 --make-bed --out cases # --make-bed could be inferred.
plink --data out.controls --oxford-single-chr 1 --make-bed --out controls
plink --bfile cases --bmerge controls --out merged # --make-bed here errors b/c of ambiguous sex.

plink --bfile merged --allow-no-sex --assoc # Normally, ambiguous-sex samples are discarded.
plink --bfile merged --allow-no-sex --logistic
