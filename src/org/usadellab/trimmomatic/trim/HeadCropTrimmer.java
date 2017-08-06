package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

public class HeadCropTrimmer extends AbstractSingleRecordTrimmer
{
	private int bases;
	private int maxLength=Integer.MAX_VALUE/2;

	public HeadCropTrimmer(String args)
	{
		String arg[]=args.split(":");
	
		bases=Integer.parseInt(arg[0]);
		
		if(arg.length>1)
			maxLength=Integer.parseInt(arg[1]);
			
	}

	@Override
	public FastqRecord processRecord(FastqRecord in)
	{
		int len=in.getSequence().length();
		
		int toTrim=bases;
		int overLen=len-toTrim-maxLength;
		
		if(overLen>0)
			toTrim+=overLen;

		if(len<=toTrim)
			return null;
		
		if(toTrim==0)
			return in;
	
		return new FastqRecord(in,toTrim,len-toTrim);
	}

}
