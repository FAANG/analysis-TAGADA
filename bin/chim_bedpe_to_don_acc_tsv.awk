
# ~sdjebali/Awk/chim_bedpe_to_don_acc_tsv.awk

# takes as input a chimeric junction file in bedpe format (5' block before 3' block) and outputs a tsv file with 4 columns
# that has two rows for each chimeric junction: one containing the 24 bp sequence surrounding the donor site and one containing
# the 24bp sequence surrounding the acceptor site, in the form chr\tstrand\t5pend\3pend, that can then be given to gem retriever
# to extract the actual 24bp sequences of those

# This may seem quite tricky since I need to look at different places according to the strand of the two blocks
################################################################################################################
# +/+
# exon GT ... AG exon
# -/-
# AC exon ... exon CT
# +/-
# exon GT ... exon CT
# -/+
# AC exon ... AG exon
# but in fact since we always need 12bp in exon and 12bp in intron it is about the same for all cases
# the only difference is that sometimes we have to go from the position -11 to the position +12 (when 
# exon is on the left) and sometimes from the position -12 to the position +11 (when exon is on the right)

BEGIN{OFS="\t"}

{
    split($7,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    chr1=a1[1]; 
    str1=a1[3]; 
    chr2=a2[1]; 
    str2=a2[3]; 
    if(a1[3]=="+")
    {
	if(a2[3]=="+")
	{
	    pos11=a1[2]-11; 
	    pos12=a1[2]+12; 
	    pos21=a2[2]-12; 
	    pos22=a2[2]+11;
	}
	else
	{
	    pos11=a1[2]-11; 
	    pos12=a1[2]+12; 
	    pos21=a2[2]-11; 
	    pos22=a2[2]+12;
	}
    }
    else
    {
	if(a2[3]=="+")
	{
	    pos11=a1[2]-12; 
	    pos12=a1[2]+11; 
	    pos21=a2[2]-12; 
	    pos22=a2[2]+11;
	}
	else
	{
	    pos11=a1[2]-12; 
	    pos12=a1[2]+11; 
	    pos21=a2[2]-11; 
	    pos22=a2[2]+12;
	}
    } 
    print chr1, str1, pos11, pos12; 
    print chr2, str2, pos21, pos22;
}