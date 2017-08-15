package org.usadellab.trimmomatic.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayDeque;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.zip.ZipInputStream;

import org.itadaki.bzip2.BZip2InputStream;
import org.usadellab.trimmomatic.util.ConcatGZIPInputStream;
import org.usadellab.trimmomatic.util.PositionTrackingInputStream;

import ngs.ErrorMsg;
import ngs.ReadCollection;
import ngs.Read;
import ngs.ReadIterator;

public class FastqParser {

	private static final int PREREAD_COUNT=1000;

    private int phredOffset;
    private ArrayDeque<FastqRecord> deque;
    int qualHistogram[];
    int patternHistogram[];
    
    private PositionTrackingInputStream posTrackInputStream;
    private BufferedReader reader1=null;
    private BufferedReader reader2=null;
    private FastqRecord[] current;
    private long fileLength;
    private long sra_max_spot;
    private long sra_processed_spot;
    private boolean single;
    private boolean sra;
    private AtomicBoolean atEOF;

    private ReadCollection run;
    private String run_name;
    private ReadIterator iter;
    
    public FastqParser(int phredOffset, File file) throws IOException {
        this.phredOffset = phredOffset;
        this.sra=false;
        this.deque=new ArrayDeque<FastqRecord>(PREREAD_COUNT);
        this.reader1=parse(file);
        this.single=true;
        this.atEOF=new AtomicBoolean();

        parseOne();

    }

    public FastqParser(int phredOffset, File  file1, File file2) throws IOException {
        this.phredOffset = phredOffset;
        this.sra=false;
        this.deque=new ArrayDeque<FastqRecord>(PREREAD_COUNT);
        this.reader1=parse(file1);
        this.reader2=parse(file2);
        this.single=false;
        this.atEOF=new AtomicBoolean();
        parseOne();
    }

    public FastqParser(int phredOffset, String acc) throws ErrorMsg, Exception {
        this.phredOffset = phredOffset;
        this.sra=true;
        run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection(acc);
        this.sra_max_spot=run.getReadCount();
        this.iter=run.getReadRange(1, this.sra_max_spot);

        this.deque=new ArrayDeque<FastqRecord>(PREREAD_COUNT);
        this.atEOF=new AtomicBoolean();
        parseOne();
    }

/*    public void setPhredOffset(int phredOffset)
    {
    	this.phredOffset=phredOffset;
    	
    	if(current!=null)
    		current.setPhredOffset(phredOffset);
    }*/
    
    public void parseOne() throws IOException 
    {
        current = null;

        String name;
        String sequence;
        String comment;
        String quality;

        if (this.sra) {
            try {
            if(iter.nextRead()){
                sra_processed_spot++;
                int n=0;
                if (iter.getNumFragments()==2) {
                    current=new FastqRecord[2];
                    single=false;
                }
                else{
                    current=new FastqRecord[1];
                    single=true;
                }
                while( iter.nextFragment() ){
                    n++;
                    //String Fid[]=it.getFragmentId().split("\\.");
                    //int i = Integer.parseInt(Fid[1].substring(2));
                    
                    name=iter.getReadName()+"/"+n;
                    sequence=iter.getFragmentBases();
                    comment=iter.getReadId();
                    quality=iter.getFragmentQualities();

                    if (n==1){
                        current[0]=new FastqRecord(name, sequence, comment, quality, this.phredOffset);
                    }
                    else {
                        current[1]=new FastqRecord(name, sequence, comment, quality, this.phredOffset) ;
                    }

                }
            }

            else atEOF.set(true);
        }
        catch (ErrorMsg x){
            System . err . println ( x . toString () );
        }
        catch (Exception x){
            System . err . println ( x . toString () );
        }
        }
        else {
            String line;
            if (single) {
                current=new FastqRecord[1];
            }
            else {
                current=new FastqRecord[2];
            }

            line = reader1.readLine();
            if (line == null) {
                atEOF.set(true);
                return;
            }
        
            if (line.charAt(0)=='@') {
                name = line.substring(1);
            } else {
                throw new RuntimeException("Invalid FASTQ name line: " + line);
            }

            sequence = reader1.readLine();
            if(sequence==null)
                throw new RuntimeException("Missing sequence line from record: " + name);

            line = reader1.readLine();
            if(line==null)
                throw new RuntimeException("Missing comment line from record: " + name);

            if (line.charAt(0)=='+') {
                comment = line.substring(1);
            } else {
                throw new RuntimeException("Invalid FASTQ comment line: " + line);
            }

            quality = reader1.readLine();
            if(quality==null)
                throw new RuntimeException("Missing quality line from record: " + name);

            current[0]=new FastqRecord(name, sequence, comment, quality, this.phredOffset);

            if (! single){
            line = reader2.readLine();
            if (line == null) {
                atEOF.set(true);
                return;
            }
        
            if (line.charAt(0)=='@') {
                name = line.substring(1);
            } else {
                throw new RuntimeException("Invalid FASTQ name line: " + line);
            }

            sequence = reader2.readLine();
            if(sequence==null)
                throw new RuntimeException("Missing sequence line from record: " + name);

            line = reader2.readLine();
            if(line==null)
                throw new RuntimeException("Missing comment line from record: " + name);

            if (line.charAt(0)=='+') {
                comment = line.substring(1);
            } else {
                throw new RuntimeException("Invalid FASTQ comment line: " + line);
            }

            quality = reader2.readLine();
            if(quality==null)
                throw new RuntimeException("Missing quality line from record: " + name);

            current[1]=new FastqRecord(name, sequence, comment, quality, this.phredOffset) ;

        }

    }
}


    public int getProgress() {
    	if(atEOF.get())
    		return 100;
    	
        if (sra){
            return (int)(((float) sra_processed_spot / sra_max_spot) * 100); 
        }
        else {
    	   long bytesRead=posTrackInputStream.getPosition();
    	
    	   return (int)(((float) bytesRead / fileLength) * 100);    
        }
    }

    
    private void accumulateHistogram(FastqRecord rec)
    {
    	int quals[]=rec.getQualityAsInteger(false);
    	
    	for(int i: quals)
    		qualHistogram[i]++;
    }
    
    public int determinePhredOffset()
    {
    	int phred33Total=0;
    	int phred64Total=0;

    	for(int i=33;i<=58;i++)
    		phred33Total+=qualHistogram[i];
    	
    	for(int i=80;i<=104;i++)
    		phred64Total+=qualHistogram[i];
    	
    	if(phred33Total==0 && phred64Total>0)
    		return 64;

    	if(phred64Total==0 && phred33Total>0)
    		return 33;
    	
    	return 0;
    }
    
    
    public BufferedReader parse(File file) throws IOException {
        String name = file.getName();
        fileLength = file.length();
        
        posTrackInputStream=new PositionTrackingInputStream(new FileInputStream(file));
        
        InputStream contentInputStream=posTrackInputStream;
        
        if (name.toLowerCase().endsWith(".gz")) {
            contentInputStream=new ConcatGZIPInputStream(posTrackInputStream);
        } else if (name.toLowerCase().endsWith(".bz2")) {
            contentInputStream=new BZip2InputStream(posTrackInputStream, false);
        } else if (name.toLowerCase().endsWith(".zip")) {
            contentInputStream=new ZipInputStream(posTrackInputStream);
        }
        
        BufferedReader reader=new BufferedReader(new InputStreamReader(contentInputStream), 32768);
        return reader;
    }

    public void close() throws IOException {
        if(reader1!=null) reader1.close();
        if(reader2!=null) reader2.close();
    }

    public boolean hasNext() {
        return (!deque.isEmpty()) || (current != null);
    }

    public FastqRecord[] next() throws IOException {

    		FastqRecord[] recs = this.current;
    		parseOne();
    		return recs;

    }

}
