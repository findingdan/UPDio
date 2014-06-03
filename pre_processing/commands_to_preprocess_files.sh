# These are the commands were used to produce the sorted and homRef-added VCF files used as input for UPDio

zcat child.exome.vcf.gz | ../scripts/sort-vcf | bgzip > child.exome.sorted.vcf.gz
zcat mom.exome.vcf.gz | ../scripts/sort-vcf | bgzip > mom.exome.sorted.vcf.gz
zcat dad.exome.vcf.gz | ../scripts/sort-vcf | bgzip > dad.exome.sorted.vcf.gz

perl ../scripts/add_hom_refs_to_vcf.pl --polymorphic_sites ../sample_data/common_variants_within_well_covered_target_regions.txt --no_homREF_vcf child.exome.sorted.vcf.gz | bgzip > child.exome.sorted.homREFed.vcf.gz
perl ../scripts/add_hom_refs_to_vcf.pl --polymorphic_sites ../sample_data/common_variants_within_well_covered_target_regions.txt --no_homREF_vcf mom.exome.sorted.vcf.gz | bgzip > mom.exome.sorted.homREFed.vcf.gz
perl ../scripts/add_hom_refs_to_vcf.pl --polymorphic_sites ../sample_data/common_variants_within_well_covered_target_regions.txt --no_homREF_vcf dad.exome.sorted.vcf.gz | bgzip > dad.exome.sorted.homREFed.vcf.gz
