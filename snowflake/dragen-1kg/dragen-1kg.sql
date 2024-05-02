-- Data location from 1000-genomes files
-- Create stage including Directory Table
create or replace stage dragen_all
   directory = (enable = true)
   url = 's3://1000genomes-dragen-3.7.6/data/individuals/hg38-graph-based'
   file_format = (type = CSV compression = AUTO)
;

-- Sample query from Directory Table showing a pattern selection of 8 Dragan files
SELECT * FROM DIRECTORY( @dragen_all )
where relative_path rlike '/HG0011[0-7].*.hard-filtered.vcf.gz'
;

-- Quick Peek (200 rows) at a single VCF file:
select * from table(ingest_vcf(BUILD_SCOPED_FILE_URL(
    @dragen_all, '/HG00116/HG00116.hard-filtered.vcf.gz'), 200));

create or replace table GENOTYPES_BY_SAMPLE (
   CHROM       varchar,
   POS         integer,
   ID          varchar,
   REF         varchar,
   ALT         array  ,
   QUAL        integer,
   FILTER      varchar,
   INFO        variant,
   SAMPLE_ID   varchar,
   VALS        variant,
   ALLELE1     varchar,
   ALLELE2     varchar,
   FILENAME    varchar
);

insert into GENOTYPES_BY_SAMPLE (
   CHROM       ,
   POS         ,
   ID          ,
   REF         ,
   ALT         ,
   QUAL        ,
   FILTER      ,
   INFO        ,
   SAMPLE_ID   ,
   VALS        ,
   ALLELE1     ,
   ALLELE2     ,
   FILENAME    )
with file_list as (
   SELECT file_url, relative_path FROM DIRECTORY( @dragen_all )
   --where relative_path rlike '.*.hard-filtered.vcf.gz'  //select all 3202 genomes
   where relative_path rlike '/HG0011[0-7].*.hard-filtered.vcf.gz'  //selection of 8 genomes
   --where relative_path rlike '/HG00[0-2][0-9][0-9].*.hard-filtered.vcf.gz'  //selection of 122 genomes
 order by random()  //temporary hint to maximize parallelism  
)
select
   replace(CHROM, 'chr','')         , 
   POS           ,
   ID            ,
   REF           ,
   split(ALT,','),
   QUAL          ,
   FILTER        ,
   INFO          , 
   SAMPLEID      ,
   SAMPLEVALS    ,
   allele1(ref, split(ALT,','), SAMPLEVALS:GT::varchar) as ALLELE1,
   allele2(ref, split(ALT,','), SAMPLEVALS:GT::varchar) as ALLELE2,
   split_part(relative_path, '/',3)
from file_list,
    table(ingest_vcf(BUILD_SCOPED_FILE_URL(@dragen_all, relative_path), 0)) vcf
;

-- Percentile distributions of Sampling Depths
select distinct sample_id,
  percentile_cont(.05) within group (order by vals:DP::number) over (partition by(sample_id)),
  percentile_cont(.1) within group (order by vals:DP::number) over (partition by(sample_id)),
  percentile_cont(.5) within group (order by vals:DP::number) over (partition by(sample_id)),
  percentile_cont(.9) within group (order by vals:DP::number) over (partition by(sample_id)),
  percentile_cont(.95) within group (order by vals:DP::number) over (partition by(sample_id))
from genotypes_by_sample
;

-- Total Count of variant calls by type, per chrom
select chrom,
   sum(case when allele1 = ref and allele2 = ref then 1 else 0 end) as homozref,
   sum(case when allele1 != allele2 then 1 else 0 end) as heterozalt,
   sum(case when allele1 != ref and allele2 = allele1 then 1 else 0 end) as homozalt
from Genotypes_by_sample
where try_to_number(chrom) is not null
group by 1 order by try_to_number(chrom);

create or replace stage clinvar
   url = 's3://aws-roda-hcls-datalake/clinvar_summary_variants'
   file_format = (type = PARQUET)
;


create or replace table CLINVAR (
 chrom                  varchar(2),
 pos                    number    ,
 ref                    varchar   ,
 alt                    varchar   ,
 ALLELEID               number    ,
 chromosomeaccession    varchar   ,        
 CLNSIG                 varchar   ,
 clinsigsimple          number    ,  
 cytogenetic            varchar   ,
 geneid                 number    ,
 genesymbol             varchar   ,
 guidelines             varchar   ,
 hgnc_id                varchar   ,
 lastevaluated          varchar   ,  
 name                   varchar   ,   
 numbersubmitters       number    ,     
 origin                 varchar   ,
 originsimple           varchar   , 
 otherids               varchar   ,
 CLNDISB                array     ,
 CLNDN                  array     ,
 rcvaccession           array     , 
 CLNREVSTAT             varchar   ,
 RS                     number    ,
 startpos               number    ,
 stoppos                number    ,
 submittercategories    number    ,        
 testedingtr            varchar   ,
 type                   varchar   ,
 variationid            number    ,
 full_annotation        variant  
);


insert into CLINVAR
select
   $1:chromosome                 ::varchar(2)  chrom,
   $1:positionvcf                ::number      pos,
   $1:referenceallelevcf         ::varchar     ref,
   $1:alternateallelevcf         ::varchar     alt,
   $1:alleleid                   ::number      ALLELEID,
   $1:chromosomeaccession        ::varchar     chromosomeaccession,
   $1:clinicalsignificance       ::varchar     CLNSIG,
   $1:clinsigsimple              ::number      clinsigsimple,
   $1:cytogenetic                ::varchar     cytogenetic,
   $1:geneid                     ::number      geneid,
   $1:genesymbol                 ::varchar     genesymbol,
   $1:guidelines                 ::varchar     guidelines,
   $1:hgnc_id                    ::varchar     hgnc_id,
   $1:lastevaluated              ::varchar     lastevaluated,
   $1:name                       ::varchar     name,   
   $1:numbersubmitters           ::number      numbersubmitters,
   $1:origin                     ::varchar     origin,
   $1:originsimple               ::varchar     originsimple,
   $1:otherids                   ::varchar     otherids,
   split($1:phenotypeids::varchar,  '|')       CLNDISB,
   split($1:phenotypelist::varchar, '|')       CLNDN,
   split($1:rcvaccession::varchar,  '|')       rcvaccession,
   $1:reviewstatus               ::varchar     CLNREVSTAT,
   $1:rsid_dbsnp                 ::number      RS,
   $1:start                      ::number      startpos,
   $1:stop                       ::number      stoppos,
   $1:submittercategories        ::number      submittercategories,
   $1:testedingtr                ::varchar     testedingtr,
   $1:type                       ::varchar     type,
   $1:variationid                ::number      variationid,
   $1                            ::variant     full_annotation
FROM '@clinvar/variant_summary/'
where $1:assembly = 'GRCh38'
order by chrom, pos
;

create or replace stage population
   url = 's3://1000genomes/1000G_2504_high_coverage/additional_698_related/'
   file_format = (type = CSV compression = AUTO field_delimiter=' ' skip_header=1)
   ;


create or replace table panel (
  Sample_ID        varchar,
  Family_ID        varchar,
  Father_ID        varchar,
  Mother_ID        varchar,
  Gender           varchar,
  Population       varchar,
  Superpopulation  varchar
);


insert into panel
select $2, $1, $3, $4, case when $5 = 1 then 'male' else 'female' end, $6, $7
from '@population/20130606_g1k_3202_samples_ped_population.txt'
;


-- Example Query:  Britain Female Samples for subset of chrom 10
select g.sample_id, chrom, pos, ref, alt, vals:GT::varchar gt,
   allele1,
   allele2
from genotypes_by_sample g
join panel p on p.sample_id = g.sample_id
where chrom = '10' and pos between 100000 and 500000
   and population = 'GBR'
   and gender = 'female'
limit 100;


-- Example Query -- locations associated with hereditary colon cancer
select   g.chrom, g.pos, g.ref, allele1, allele2, count (sample_id)
from genotypes_by_sample g
join clinvar c
   on c.chrom = g.chrom and c.pos = g.pos and c.ref = g.ref
where
   array_contains('Hereditary nonpolyposis colorectal neoplasms'::variant, CLNDN)
group by 1,2,3,4,5
order by 1,2,5
;


-- Example Query -- specific variants associated with hereditary colon cancer
select   g.chrom, g.pos, g.ref, c.alt clinvar_alt, genesymbol, allele1, allele2, count (sample_id)
from genotypes_by_sample g
join clinvar c
   on c.chrom = g.chrom and c.pos = g.pos and c.ref = g.ref
   and ((Allele1= c.alt) or (Allele2= c.alt))
where
   array_contains('Hereditary nonpolyposis colorectal neoplasms'::variant, CLNDN)
group by 1,2,3,4,5,6,7
order by 1,2,5
;