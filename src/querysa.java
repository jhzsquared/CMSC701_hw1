import java.io.*;
import java.util.*;
import java.util.Collections;

public class querysa {
    // read in binary string and suffix array (and maybe prefix table) 
    static SuffixArrayProto.SuffixMessage ReadSABin(String binPath){
        // Read in protobuf output
        SuffixArrayProto.SuffixMessage mess=null;
        try (FileInputStream input = new FileInputStream(binPath)) {
            SuffixArrayProto.SuffixMessage.Builder messageRead = SuffixArrayProto.SuffixMessage.newBuilder();
            mess = messageRead.mergeFrom(input).build();

        } catch (IOException e) {
            System.out.println(e);
            System.exit(1);
        }
        return mess;
    }

    static Map<String,String> ReadFastaQuery(String queryPath){
         //flie path needs to be full path   
         // assumes onee sequence per file 
        String query,header = "";
        boolean firstLine = true;
        StringBuilder stringBuilder = new StringBuilder();
        Map<String,String> queries = new HashMap<String, String>();

         try (BufferedReader reader = new BufferedReader(new FileReader(queryPath))){
             String line = reader.readLine();
             while (line != null) {
                 if (line.charAt(0) == '>' && firstLine) {
                    //new query coming up
                     stringBuilder = new StringBuilder();
                     String[] splitHeader = line.split("\\s+");
                     header = splitHeader[0].substring(1); // remove the >
                     firstLine = false;
                 } else if (line.charAt(0) != '>'){
                    //getting query lines
                     stringBuilder.append(line.toUpperCase());
                 } else if (line.charAt(0) == '>' && firstLine==false){
                    //new query request, so save previous reesults
                    query = stringBuilder.toString();
                    queries.put(header,query);
                    //initiate for next query
                    stringBuilder = new StringBuilder();
                    String[] splitHeader = line.split("\\s+");
                    header = splitHeader[0].substring(1); // remove the >
                 }
                 line = reader.readLine();
             }
             reader.close();
             //save final query result
            query = stringBuilder.toString();
            queries.put(header,query);
 
         } catch (FileNotFoundException e){
             System.err.println(e);
             System.exit(1);
          } catch (IOException e)  {
             System.out.println(e);
             System.exit(1);
         }
         return queries;
    }

    static String getSubstring(String genome, int start, int length){
        String substring;
        try{
            substring = genome.substring(start, start + length);
        } catch (java.lang.StringIndexOutOfBoundsException e){ //start is too short
            substring = genome.substring(start);
        }
        return substring;
    }
    
    static List<Integer> naiveQuery(String query, String genome, List<Integer> suffixArray, Map<String,SuffixArrayProto.indexInterval> prefixTable){
        // naive binary search (if prefix table is not null then it will use it)
        List<Integer> hits = new ArrayList<Integer>();
        int l= 0;
        int r = suffixArray.size();
        int queryLength = query.length();
        int c = (int)Math.floor((l+r)/2);
        //start timing
        long startTime = System.nanoTime();

        if (prefixTable.size()!=0){
        // set initial left and right limits based on prefixtable findings
            //check if pattern is longer or shorter than prefixes
            int prefixLength = prefixTable.keySet().iterator().next().length();
            if (prefixLength > queryLength){
                //query is shorter than prefixes generated
                //find range of prefixes that fit //binary search --and then you have effectively found range
                //binary search for lower bound
                List<String> prefixList = new ArrayList<String>(prefixTable.keySet());
                Collections.sort(prefixList); //lexicographically sort
                l = 0;
                r = prefixList.size();
                while (l< r){
                    c = (int)Math.floor((l+r)/2);
                    //check middle
                    if (getSubstring(prefixList.get(c), 0, queryLength).compareTo(query) >= 0){
                        // query is  smaller or equal  to center substring
                        r = c; //left bisect
                    } else{
                        //query is bigger than center
                        l = c+1;// right bisect
                    }
                }
                int leftBound = l;
                //binary search for upper bound (might as well start at the left bound)
                r = prefixList.size();
                while (l<r){
                    c = (int)Math.floor((l+r)/2);
                    //check middle
                    if (getSubstring(prefixList.get(c), 0, queryLength).compareTo(query) <= 0){
                        // query is bigger or equal  to center substring
                        l = c+1; //right bisect
                    } else{
                        //query is bigger than center
                        r = c;// right bisect
                    }
                }
                // convert prefixtable hits to actual genome indices
                leftBound = prefixTable.get(prefixList.get(leftBound)).getIndices(0);
                r = prefixTable.get(prefixList.get(r)).getIndices(1);
                for (int i = leftBound; i<r; i++){
                    hits.add(suffixArray.get(i).intValue());
                }   
                //stop timing
                long endTime = System.nanoTime();
                long duration = (endTime-startTime);
                System.out.println("Query of "+ query.length() + " length took " + duration);
                return hits;
                
            } else { 
                //query is longer than prefix, proceed with standard search
                String queryPrefix = query.substring(0, prefixLength);
                if (prefixTable.get(queryPrefix)==null){
                    return hits; //no matches if prefix can not be found
                } else{
                    l = prefixTable.get(queryPrefix).getIndices(0);
                    r = prefixTable.get(queryPrefix).getIndices(1);
                }             
            }
        }
        //binary search for lower bound
        while (l< r){
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue(),queryLength).compareTo(query) >= 0){
                // query is  smaller or equal  to center substring
                r = c; //left bisect
            } else{
                //query is bigger than center
                l = c+1;// right bisect
            }
        }
        int leftBound = l;
        //binary search for upper bound (might as well start at the left bound)
        r = suffixArray.size();
        while (l<r){
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue(),queryLength).compareTo(query) <= 0){
                // query is bigger or equal  to center substring
                l = c+1; //right bisect
            } else{
                //query is bigger than center
                r = c;// right bisect
            }
        }
        // convert suffix Array hits to actual genome indices
        for (int i = leftBound; i<r; i++){
            hits.add(suffixArray.get(i).intValue());
        }   
        //stop timing
        long endTime = System.nanoTime();
        long duration = (endTime-startTime);
        System.out.println("Query of "+ query.length() + " length took " + duration);
        return hits;
    }
    static List<Integer> simpAccelQuery(String query, String genome, List<Integer> suffixArray, Map<String,SuffixArrayProto.indexInterval> prefixTable){
        List<Integer> hits = new ArrayList<Integer>();
        int l= 0;
        int r = suffixArray.size();
        int queryLength = query.length();
        int c = (int)Math.floor((l+r)/2);
        //start timing
        long startTime = System.nanoTime();

        if (prefixTable.size()!=0){
        // set initial left and right limits based on prefixtable findings
            //check if pattern is longer or shorter than prefixes
            int prefixLength = prefixTable.keySet().iterator().next().length();
            if (prefixLength > queryLength){
                //query is shorter than prefixes generated
                //find range of prefixes that fit //binary search --and then you have effectively found range
                //binary search for lower bound
                List<String> prefixList = new ArrayList<String>(prefixTable.keySet());
                Collections.sort(prefixList); //lexicographically sort
                l = 0;
                r = prefixList.size();
                while (l< r){
                    c = (int)Math.floor((l+r)/2);
                    //check middle
                    if (getSubstring(prefixList.get(c), 0, queryLength).compareTo(query) >= 0){
                        // query is  smaller or equal  to center substring
                        r = c; //left bisect
                    } else{
                        //query is bigger than center
                        l = c+1;// right bisect
                    }
                }
                int leftBound = l;
                //binary search for upper bound (might as well start at the left bound)
                r = prefixList.size();
                while (l<r){
                    c = (int)Math.floor((l+r)/2);
                    //check middle
                    if (getSubstring(prefixList.get(c), 0, queryLength).compareTo(query) <= 0){
                        // query is bigger or equal  to center substring
                        l = c+1; //right bisect
                    } else{
                        //query is bigger than center
                        r = c;// right bisect
                    }
                }
                // convert prefixtable hits to actual genome indices
                leftBound = prefixTable.get(prefixList.get(leftBound)).getIndices(0);
                r = prefixTable.get(prefixList.get(r)).getIndices(1);
                for (int i = leftBound; i<r; i++){
                    hits.add(suffixArray.get(i).intValue());
                }   
                //stop timing
                long endTime = System.nanoTime();
                long duration = (endTime-startTime);
                System.out.println("Query of "+ query.length() + " length took " + duration);
                return hits;
                
            } else { 
                //query is longer than prefix, proceed with standard search
                String queryPrefix = query.substring(0, prefixLength);
                if (prefixTable.get(queryPrefix)==null){
                    return hits; //no matches if prefix can not be found
                } else{
                    l = prefixTable.get(queryPrefix).getIndices(0);
                    r = prefixTable.get(queryPrefix).getIndices(1);
                }             
            }
        }
        //binary search for lower bound
        while (l< r){
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue(),queryLength).compareTo(query) >= 0){
                // query is  smaller or equal  to center substring
                r = c; //left bisect
            } else{
                //query is bigger than center
                l = c+1;// right bisect
            }
        }
        int leftBound = l;
        //binary search for upper bound (might as well start at the left bound)
        r = suffixArray.size();
        while (l<r){
            c = (int)Math.floor((l+r)/2);
            //check middle
            if (getSubstring(genome, suffixArray.get(c).intValue(),queryLength).compareTo(query) <= 0){
                // query is bigger or equal  to center substring
                l = c+1; //right bisect
            } else{
                //query is bigger than center
                r = c;// right bisect
            }
        }
        // convert suffix Array hits to actual genome indices
        for (int i = leftBound; i<r; i++){
            hits.add(suffixArray.get(i).intValue());
        }   
        //stop timing
        long endTime = System.nanoTime();
        long duration = (endTime-startTime);
        System.out.println("Query of "+ query.length() + " length took " + duration);
        return hits;
    }

    public static void main(String[] args) {
        String binPath, queryPath, queryMode, output;
        boolean preftab;
        if (args.length==4){
            binPath = args[0]; //path to binary file
            queryPath = args[1]; //path to query fasta file
            queryMode = args[2]; //naive or simpaccel
            output = args[3]; //file path for output
        } else {
            System.err.println("Invalid arguments count: " + args.length);
            return;
        }
        // read in protobuf suffix array
        SuffixArrayProto.SuffixMessage suffixMessage = ReadSABin(binPath);
        String genome = suffixMessage.getGenome();
        System.out.println("Reference is " + genome.length() + " characters long");
        List<Integer> suffixArray = suffixMessage.getSuffixArrayList(); //suffixArray.get(i).intValue()
        Map<String,SuffixArrayProto.indexInterval> prefixTable = suffixMessage.getPreftabMap(); //prefixTable.get('pattern').getIndices(index)
        if (prefixTable.size()==0){
            preftab = false;
        } else {
            preftab = true;
        }
       
        // read in fasta query
        Map<String,String> queries = ReadFastaQuery(queryPath);
        //conduct query
        for (Map.Entry<String, String> entry : queries.entrySet()){
            List<Integer> queryResults = naiveQuery(entry.getValue(), genome, suffixArray, prefixTable);
            System.out.println(queryResults);
        }
       

        // save output
    }
}
