vt decompose -s                         a.vcf -o a-vt.vcf
bcftools norm -m-any                    a.vcf -o a-bcftools_m.vcf
bcftools norm -m-any --multi-overlaps . a.vcf -o a-bcftools_m_dot.vcf
bcftools norm -a                        a.vcf -o a-bcftools_a.vcf
bcftools norm -a --atom-overlaps '*'    a.vcf -o a-bcftools_a_star.vcf
bcftools norm -a --atom-overlaps .      a.vcf -o a-bcftools_a_dot.vcf

vt decompose -s                         b.vcf -o b-vt.vcf
