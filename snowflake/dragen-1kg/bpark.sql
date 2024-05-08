SELECT concat('aws s3 cp --no-sign-request s3://1000genomes-dragen-3.7.6/data/individuals/hg38-graph-based/', relative_path, ' ./1kg/') AS cli FROM DIRECTORY( @dragen_all )
where relative_path rlike '/*.*.hard-filtered.vcf.gz'
ORDER BY SIZE DESC 
LIMIT 1000 ;