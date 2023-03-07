package textproc;

import java.util.Map;
import java.util.TreeMap;


public class MultiWordCounter implements TextProcessor {
	Map<String, Integer> m;
	
	public MultiWordCounter(String [] s) {
		m = new TreeMap<String, Integer>();
		for (int i = 0; i < s.length; i++) {
			m.put(s[i], 0);
		}
	}
	
	public void process(String w) {
		if(m.containsKey(w)) {
			m.replace(w, m.get(w) + 1);
		}
	}
	
	public void report() {
		for (String key : m.keySet()) {
			System.out.println(key + ": " + m.get(key));
		}
	}
	
}
