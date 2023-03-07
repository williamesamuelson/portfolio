package textproc;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class GeneralWordCounter implements TextProcessor {
	private Map<String, Integer> m;
	private Set<String> stopwords;

	public GeneralWordCounter(Set<String> stopwords) {
		m = new TreeMap<String, Integer>();
		this.stopwords = stopwords;
	}

	public void process(String w) {
		if (!stopwords.contains(w)) {
			if (m.containsKey(w)) {
				m.replace(w, m.get(w) + 1);
			} else {
				m.put(w, 1);
			}
		}
	}

	public void report() {
//		for (String key : m.keySet()) {
//			if (m.get(key) >= 200) {
//				System.out.println(key + ": " + m.get(key));
//			}
//		}

		Set<Map.Entry<String, Integer>> wordSet = m.entrySet();
		List<Map.Entry<String, Integer>> wordList = new ArrayList<>(wordSet);
		wordList.sort(new WordCountComparator());

		for (int i = 0; i < 20; i++) {
			System.out.println(wordList.get(i).getKey() + ": " + wordList.get(i).getValue());
		}
	}

	public Set<Map.Entry<String, Integer>> getWords() {
		return m.entrySet();
	}
}
