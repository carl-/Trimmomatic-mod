package org.usadellab.trimmomatic;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.fastq.FastqParser;
import org.usadellab.trimmomatic.fastq.FastqRecord;
//import org.usadellab.trimmomatic.fastq.FastqRecordWrapper;
import org.usadellab.trimmomatic.fastq.FastqSerializer;
import org.usadellab.trimmomatic.fastq.PairingValidator;
import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ParserWorker;
import org.usadellab.trimmomatic.threading.SerializerWorker;
import org.usadellab.trimmomatic.threading.TrimLogWorker;
import org.usadellab.trimmomatic.threading.TrimStatsWorker;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.trim.TrimmerFactory;
import org.usadellab.trimmomatic.util.Logger;

import ngs.ErrorMsg;

public class TrimmomaticPE extends Trimmomatic
{

	/**
	 * Trimmomatic: The FASTQ trimmer
	 * 
	 * CROP:<LENGTH> Crop the read to specified length, by cutting off the right
	 * end LEADING:<QUAL> Trim the read by cutting off the left end while below
	 * specified quality TRAILING:<QUAL> Trim the read by cutting off the right
	 * end while below specified quality SLIDINGWINDOW:<QUAL>:<COUNT> Trim the
	 * read once the total quality of COUNT bases drops below QUAL, then trim
	 * trailing bases below QUAL
	 * 
	 * MINLEN:<LENGTH> Drop the read if less than specified length
	 */

	private Logger logger;

	public TrimmomaticPE(Logger logger)
	{	
		this.logger=logger;
	}

	public void processSingleThreaded(FastqParser parser, FastqSerializer serializer1P,
			FastqSerializer serializer1U, FastqSerializer serializer2P, FastqSerializer serializer2U,
			Trimmer trimmers[], PrintStream trimLogStream, PairingValidator pairingValidator) throws IOException
	{
		TrimStats stats = new TrimStats();

		FastqRecord originalRecs[] = new FastqRecord[2];
		FastqRecord recs[] = new FastqRecord[2];

		while (parser.hasNext())
			{
			FastqRecord[] r;
			//try {
				r= parser.next();
					if(recs.length!=2){
					System.err.println(" Expect Paired Reads, but get only one Read !");
					throw new RuntimeException("Invalid Paired Reads");
				}
				originalRecs[0] = recs[0] = r[0];
				originalRecs[1] = recs[1] = r[1];
			//}
			//catch ( ErrorMsg x){
			//	System . err . println ( x . toString () );
			//}
			//catch (Exception x){
			//	System . err . println ( x . toString () );
			//	throw x;
			//}


			if(pairingValidator!=null)
				pairingValidator.validatePair(recs[0], recs[1]);
			
			for (int i = 0; i < trimmers.length; i++)
				{
				try
					{
					recs = trimmers[i].processRecords(recs);
					}
				catch (RuntimeException e)
					{
					logger.errorln("Exception processing reads: " + originalRecs[0].getName() + " and "
							+ originalRecs[1].getName());
					e.printStackTrace();
					throw e;
					}
				}

			if (recs[0] != null && recs[1] != null)
				{
				serializer1P.writeRecord(recs[0]);
				serializer2P.writeRecord(recs[1]);
				}
			else if (recs[0] != null)
				serializer1U.writeRecord(recs[0]);
			else if (recs[1] != null)
				serializer2U.writeRecord(recs[1]);

			stats.logPair(originalRecs, recs);

			if (trimLogStream != null)
				{
				for (int i = 0; i < originalRecs.length; i++)
					{
					int length = 0;
					int startPos = 0;
					int endPos = 0;
					int trimTail = 0;

					if (recs[i] != null)
						{
						length = recs[i].getSequence().length();
						startPos = recs[i].getHeadPos();
						endPos = length + startPos;
						trimTail = originalRecs[i].getSequence().length() - endPos;
						}

					trimLogStream.printf("%s %d %d %d %d\n", originalRecs[i].getName(), length, startPos, endPos,
							trimTail);
					}
				}
			}

		logger.infoln(stats.getStatsPE());
	}

	public void processMultiThreaded(FastqParser parser, FastqSerializer serializer1P,
			FastqSerializer serializer1U, FastqSerializer serializer2P, FastqSerializer serializer2U,
			Trimmer trimmers[], PrintStream trimLogStream, PairingValidator pairingValidator, int threads) throws IOException
	{
		ArrayBlockingQueue<Runnable> taskQueue = new ArrayBlockingQueue<Runnable>(threads);

		ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue1P = new ArrayBlockingQueue<Future<BlockOfRecords>>(
				threads);
		ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue1U = new ArrayBlockingQueue<Future<BlockOfRecords>>(
				threads);
		ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue2P = new ArrayBlockingQueue<Future<BlockOfRecords>>(
				threads);
		ArrayBlockingQueue<Future<BlockOfRecords>> serializerQueue2U = new ArrayBlockingQueue<Future<BlockOfRecords>>(
				threads);

		ThreadPoolExecutor taskExec = new ThreadPoolExecutor(threads, threads, 0, TimeUnit.SECONDS, taskQueue);

		SerializerWorker serializerWorker1P = new SerializerWorker(serializer1P, serializerQueue1P, 0);
		SerializerWorker serializerWorker1U = new SerializerWorker(serializer1U, serializerQueue1U, 1);
		SerializerWorker serializerWorker2P = new SerializerWorker(serializer2P, serializerQueue2P, 2);
		SerializerWorker serializerWorker2U = new SerializerWorker(serializer2U, serializerQueue2U, 3);

		Thread serializer1PThread = new Thread(serializerWorker1P);
		Thread serializer1UThread = new Thread(serializerWorker1U);
		Thread serializer2PThread = new Thread(serializerWorker2P);
		Thread serializer2UThread = new Thread(serializerWorker2U);

		ArrayBlockingQueue<Future<BlockOfRecords>> trimStatsQueue = new ArrayBlockingQueue<Future<BlockOfRecords>>(threads * 5);
		TrimStatsWorker statsWorker = new TrimStatsWorker(trimStatsQueue);
		Thread statsThread = new Thread(statsWorker);

		ArrayBlockingQueue<Future<BlockOfRecords>> trimLogQueue = null;
		TrimLogWorker trimLogWorker = null;
		Thread trimLogThread = null;

		if (trimLogStream != null)
			{
			trimLogQueue = new ArrayBlockingQueue<Future<BlockOfRecords>>(threads * 5);
			trimLogWorker = new TrimLogWorker(trimLogStream, trimLogQueue);
			trimLogThread = new Thread(trimLogWorker);
			trimLogThread.start();
			}

		serializer1PThread.start();
		serializer1UThread.start();
		serializer2PThread.start();
		serializer2UThread.start();

		statsThread.start();

		boolean done1 = false, done2 = false;

		List<FastqRecord> recs1 = null;
		List<FastqRecord> recs2 = null;
		
		try
		{
			while(parser.hasNext())
			{
			FastqRecord[] recs=parser.next();
			if(recs.length!=2){
				System.err.println(" Expect Paired Reads, but get only one Read !");
				throw new RuntimeException("Invalid Paired Reads");
			}
			recs1=new ArrayList<FastqRecord>();
			recs2=new ArrayList<FastqRecord>();
			recs1.add(recs[0]);
			recs2.add(recs[1]);
			if(pairingValidator!=null)
				pairingValidator.validatePairs(recs1, recs2);
			
			BlockOfRecords bor = new BlockOfRecords(recs1,recs2);
			BlockOfWork work = new BlockOfWork(logger, trimmers, bor, true, trimLogStream != null);
			while (taskQueue.remainingCapacity() < 1)
				Thread.sleep(100);

			Future<BlockOfRecords> future = taskExec.submit(work);
			
			serializerQueue1P.put(future);
			serializerQueue1U.put(future);
			serializerQueue2P.put(future);
			serializerQueue2U.put(future);
			trimStatsQueue.put(future);

			if (trimLogQueue != null)
				trimLogQueue.put(future);
			}
			BlockOfWork work = new BlockOfWork(logger, trimmers, new BlockOfRecords(new ArrayList<FastqRecord>(),new ArrayList<FastqRecord>()), true, trimLogStream != null);
			Future<BlockOfRecords> future = taskExec.submit(work);
			serializerQueue1P.put(future);
			serializerQueue1U.put(future);
			serializerQueue2P.put(future);
			serializerQueue2U.put(future);
			trimStatsQueue.put(future);
			if (trimLogQueue != null)
				trimLogQueue.put(future);
		
			parser.close();
		
			taskExec.shutdown();
			taskExec.awaitTermination(1, TimeUnit.HOURS);
		
			serializer1PThread.join();
			serializer1UThread.join();
			serializer2PThread.join();
			serializer2UThread.join();

			if (trimLogThread != null){
				trimLogThread.join();
				statsThread.join();
				logger.infoln(statsWorker.getStats().getStatsPE());
			}
		}
		catch (InterruptedException e)
			{
			throw new RuntimeException(e);
			}
		catch (Exception e){
			System . err . println ( e . toString () );
		}
	}

	public void process(File input1, File input2, File output1P, File output1U, File output2P, File output2U,
			Trimmer trimmers[], int phredOffset, File trimLog, boolean validatePairing, int threads) throws IOException
	{
		int phred1=33;
		if(phredOffset==0)
			{
			logger.infoln("Quality encoding set to "+phred1);
			}
		else {
			phred1=phredOffset;
		}	
		FastqParser parser = new FastqParser(phred1, input1, input2);
		
		FastqSerializer serializer1P = new FastqSerializer();
		serializer1P.open(output1P);

		FastqSerializer serializer1U = new FastqSerializer();
		serializer1U.open(output1U);

		FastqSerializer serializer2P = new FastqSerializer();
		serializer2P.open(output2P);

		FastqSerializer serializer2U = new FastqSerializer();
		serializer2U.open(output2U);

		PrintStream trimLogStream = null;
		if (trimLog != null)
			// trimLogStream=new PrintStream(new BufferedOutputStream(new
			// FileOutputStream(trimLog),1000000),false);
			trimLogStream = new PrintStream(trimLog);

		PairingValidator pairingValidator=null;
		
		if(validatePairing)
			pairingValidator=new PairingValidator(logger);
		
		if (threads == 1)
			processSingleThreaded(parser, serializer1P, serializer1U, serializer2P, serializer2U, trimmers,
					trimLogStream, pairingValidator);
		else
			processMultiThreaded(parser,  serializer1P, serializer1U, serializer2P, serializer2U, trimmers,
					trimLogStream, pairingValidator, threads);

		serializer1P.close();
		serializer1U.close();
		serializer2P.close();
		serializer2U.close();

		if (trimLogStream != null)
			trimLogStream.close();
	}
	

public void process(String acc, File output1P, File output1U, File output2P, File output2U,
			Trimmer trimmers[], int phredOffset, File trimLog, boolean validatePairing, int threads) throws Exception, ErrorMsg
	{
		
		int phred1=33;
		if(phredOffset==0)
			{
			logger.infoln("Quality encoding set to "+phred1);
			}
		else {
			phred1=phredOffset;
		}	
		FastqParser parser = new FastqParser(phred1, acc);

		FastqSerializer serializer1P = new FastqSerializer();
		serializer1P.open(output1P);

		FastqSerializer serializer1U = new FastqSerializer();
		serializer1U.open(output1U);

		FastqSerializer serializer2P = new FastqSerializer();
		serializer2P.open(output2P);

		FastqSerializer serializer2U = new FastqSerializer();
		serializer2U.open(output2U);

		PrintStream trimLogStream = null;
		if (trimLog != null)
			// trimLogStream=new PrintStream(new BufferedOutputStream(new
			// FileOutputStream(trimLog),1000000),false);
			trimLogStream = new PrintStream(trimLog);

		PairingValidator pairingValidator=null;
		
		if(validatePairing)
			pairingValidator=new PairingValidator(logger);
		
		if (threads == 1)
			processSingleThreaded(parser, serializer1P, serializer1U, serializer2P, serializer2U, trimmers,
					trimLogStream, pairingValidator);
		else
			processMultiThreaded(parser, serializer1P, serializer1U, serializer2P, serializer2U, trimmers,
					trimLogStream, pairingValidator, threads);

		serializer1P.close();
		serializer1U.close();
		serializer2P.close();
		serializer2U.close();

		if (trimLogStream != null)
			trimLogStream.close();
	}
	
	private static int getFileExtensionIndex(String str)
	{
		String extensions[]={".fq",".fastq",".txt",".gz",".bz2",".zip"};
	
		String tmp=str;
		boolean done=false;
		
		while(!done)
			{
			done=true;
			for(String ext: extensions)
				{		
				if(tmp.endsWith(ext))
					{
					tmp=tmp.substring(0,tmp.length()-ext.length());
					done=false;
					}
				}
			}
	
		return tmp.length();
	}



	private static String replaceLast(String str, String out, String in)
	{
		int idx1=str.lastIndexOf(out);
		if(idx1==-1)
			return null;
		
		int idx2=idx1+out.length();
		
		return str.substring(0,idx1)+in+str.substring(idx2);
	}
	
	
	
	private static File[] calculateTemplatedInput(String baseStr)
	{
		String translation[][]={{"_R1_","_R2_"},{"_f","_r"},{".f",".r"},{"_1","_2"},{".1",".2"}};
	
		File fileBase=new File(baseStr);
		File baseDir=fileBase.getParentFile();
		
		String baseName=fileBase.getName();
		int extSplit=getFileExtensionIndex(baseName);
		
		String core=baseName.substring(0,extSplit);
		String exts=baseName.substring(extSplit);
		
		for(String pair[]: translation)
			{
			String tmp=replaceLast(core, pair[0], pair[1]);
			if(tmp!=null)
				return new File[] {fileBase, new File(baseDir, tmp+exts)};
			}
		
		return null;
	}

	
	private static File[] calculateTemplatedOutput(String baseStr)
	{
		File fileBase=new File(baseStr);
		File baseDir=fileBase.getParentFile();
		
		String baseName=fileBase.getName();
		int extSplit=getFileExtensionIndex(baseName);
		
		String core=baseName.substring(0,extSplit);
		String exts=baseName.substring(extSplit);
		
		return new File[] {new File(baseDir,core+"_1P"+exts),new File(baseDir,core+"_1U"+exts),new File(baseDir,core+"_2P"+exts),new File(baseDir,core+"_2U"+exts)};		
	}

	

	public static boolean run(String[] args) throws IOException
	{
		int argIndex = 0;
		int phredOffset = 0;
		int threads = 0;

		String templateInput=null;
		String templateOutput=null;
		
		boolean badOption = false;
		boolean validatePairs = false;
		boolean quiet=false;
		boolean showVersion=false;
		String sra=null;
		
		File trimLog = null;

		List<String> nonOptionArgs=new ArrayList<String>();
		
		while (argIndex < args.length)
			{
			String arg = args[argIndex++];
			
			if(arg.startsWith("-"))
				{			
				if (arg.equals("-phred33"))
					phredOffset = 33;
				else if (arg.equals("-phred64"))
					phredOffset = 64;
				else if (arg.equals("-threads"))
					threads = Integer.parseInt(args[argIndex++]);
				else if (arg.equals("-trimlog"))
					{
					if (argIndex < args.length)
						trimLog = new File(args[argIndex++]);
					else
						badOption = true;
					}
				else if (arg.equals("-basein"))
					{
					if (argIndex < args.length)
						templateInput = args[argIndex++];  
					else
						badOption = true;
					}
				else if (arg.equals("-baseout"))
					{
					if (argIndex < args.length)
						templateOutput = args[argIndex++];
					else
						badOption = true;
					}
				else if (arg.equals("-validatePairs"))
					validatePairs=true;
				else if (arg.equals("-quiet"))
					quiet=true;
				else if (arg.equals("-version"))
					showVersion=true; 
				else if (arg.equals("-sra"))
					{
					if (argIndex < args.length)
						sra = args[argIndex++];
					else
						badOption = true;
					}	
				else
					{
					System.err.println("Unknown option " + arg);
					badOption = true;
					}
				}
			else
				nonOptionArgs.add(arg);
			}

		if(showVersion)
			Trimmomatic.showVersion();
		
		
		int additionalArgs= sra==null ? 1+(templateInput==null?2:0)+(templateOutput==null?4:0) : 1+(templateOutput==null?4:0);
		
		if ((nonOptionArgs.size() < additionalArgs) || badOption)
			{System.err.println(nonOptionArgs.size());System.err.println(additionalArgs);return showVersion; }
		
		Logger logger=new Logger(true,true,!quiet);
		
		if(phredOffset==0) phredOffset=33;	
		logger.infoln("TrimmomaticPE: Started with arguments:");
		for (String arg : args)
			logger.info(" " + arg);
		logger.infoln();
		
		if(threads==0)
			{
			threads=calcAutoThreadCount();
			if(threads>1)
				logger.infoln("Multiple cores found: Using "+threads+" threads");
			}

		
		Iterator<String> nonOptionArgsIter=nonOptionArgs.iterator();

		File inputs[]=new File[4];
		File outputs[]=new File[4];
		if(sra==null)
		{
		if(templateInput!=null )
			{
			inputs=calculateTemplatedInput(templateInput);
			if(inputs==null)
				{
				logger.errorln("Unable to determine input files from: "+templateInput);
				System.exit(1);
				}
			
			logger.infoln("Using templated Input files: "+inputs[0]+" "+inputs[1]);
			}
		else if(templateInput==null)
			{
			inputs=new File[2];
			inputs[0]=new File(nonOptionArgsIter.next());
			inputs[1]=new File(nonOptionArgsIter.next());
			}
		
		if(templateOutput!=null)
			{
			outputs=calculateTemplatedOutput(templateOutput);
			if(outputs==null)
				{
				System.err.println("Unable to determine output files from: "+templateInput);
				System.exit(1);
				}
			
			logger.infoln("Using templated Output files: "+outputs[0]+" "+outputs[1]+" "+outputs[2]+" "+outputs[3]);
			}
		else
			{
			outputs=new File[4];
			outputs[0]=new File(nonOptionArgsIter.next());
			outputs[1]=new File(nonOptionArgsIter.next());
			outputs[2]=new File(nonOptionArgsIter.next());
			outputs[3]=new File(nonOptionArgsIter.next());
			}

		Trimmer trimmers[]=createTrimmers(logger, nonOptionArgsIter);
				
		TrimmomaticPE tm = new TrimmomaticPE(logger);

			tm.process(inputs[0], inputs[1], outputs[0], outputs[1], outputs[2], outputs[3], trimmers, phredOffset, trimLog, validatePairs, threads);
		}
		else {
					if(templateOutput!=null)
			{
			outputs=calculateTemplatedOutput(templateOutput);
			if(outputs==null)
				{
				System.err.println("Unable to determine output files from: "+templateInput);
				System.exit(1);
				}
			
			logger.infoln("Using templated Output files: "+outputs[0]+" "+outputs[1]+" "+outputs[2]+" "+outputs[3]);
			}
		else
			{
			outputs=new File[4];
			outputs[0]=new File(nonOptionArgsIter.next());
			outputs[1]=new File(nonOptionArgsIter.next());
			outputs[2]=new File(nonOptionArgsIter.next());
			outputs[3]=new File(nonOptionArgsIter.next());
			}
			Trimmer trimmers[]=createTrimmers(logger, nonOptionArgsIter);
				
			TrimmomaticPE tm = new TrimmomaticPE(logger);
			try {
				tm.process(sra, outputs[0], outputs[1], outputs[2], outputs[3], trimmers, phredOffset, trimLog, validatePairs, threads);
			}
			catch ( ErrorMsg x) {
				System . err . println ( x . toString () );
			}
			catch (Exception x){
				System . err . println ( x . toString () );
				System.exit(1);
			}
		}

		logger.infoln("TrimmomaticPE: Completed successfully");
		return true;
	}

	public static void main(String[] args) throws IOException
	{
		if (!run(args))
			{
			System.err
					.println("Usage: TrimmomaticPE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...");
			System.exit(1);
			}
	}

}
