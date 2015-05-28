#!/usr/bin/perl

#perl vcfANDMAF.pl /storage/ref/var/00-All.vcf.MAF VCF.FILES

$DBSNP = $ARGV[0];
$FILE = $ARGV[1];


open (dbSNP, $DBSNP);
open (vcfFILE, $FILE);

@dbsnp = <dbSNP>;
@vcf = <vcfFILE>;

foreach $db (@dbsnp)
{
	chomp ($db);
	@tmp = split ("\t", $db);
	$dbsnp_MAF {$tmp[0]} = $tmp[1];
}

foreach $line_vcf (@vcf)
{
	chomp ($line_vcf);
	@tmp2 = split ("\t", $line_vcf);
	
	if ($dbsnp_MAF {$tmp2[2]})
	{
		print $line_vcf."\t".$dbsnp_MAF {$tmp2[2]}."\n";
	}	
	else
	{
		print $line_vcf."\t"."."."\n";
	} 
}

close dbSNP;
close vcfFILE;
exit;
