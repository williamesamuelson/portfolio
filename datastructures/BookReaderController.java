package lab3;


import java.io.File;
import java.util.HashSet;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import javafx.application.Application;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.ListView;
import javafx.scene.control.SelectionModel;
import javafx.scene.control.TextField;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.stage.Stage;
import textproc.GeneralWordCounter;


public class BookReaderController extends Application {
	
	@Override
	public void start(Stage primaryStage) throws Exception {
		
		Scanner scan = new Scanner(new File("undantagsord.txt"));
		Set<String> stopwords = new HashSet<String>();
		while(scan.hasNext()) {
			stopwords.add(scan.next().toLowerCase());
		}
		scan.close();

		GeneralWordCounter t = new GeneralWordCounter(stopwords);
		
		Scanner s = new Scanner(new File("nilsholg.txt"));
		s.findWithinHorizon("\uFEFF", 1);
		s.useDelimiter("(\\s|,|\\.|:|;|!|\\?|'|\\\")+");
		while (s.hasNext()) {
			String word = s.next().toLowerCase();
			t.process(word);
		}
		s.close();
		
		ObservableList<Map.Entry<String, Integer>> words = FXCollections.observableArrayList(t.getWords());
		ListView<Map.Entry<String, Integer>> listView = new ListView<Map.Entry<String, Integer>>(words);
		
		BorderPane root = new BorderPane();
		Scene scene = new Scene(root, 500, 500);
		root.setCenter(listView);
		
		HBox hbox = new HBox();
		Button b1 = new Button("Alphabetic");
		Button b2 = new Button("Frequency");
		Button b3 = new Button("Find");
		b3.setDefaultButton(true);
		//Button b4 = new Button("New File...");
		TextField tf = new TextField();
	
		b1.setOnAction(event -> {
			words.sort((Map.Entry<String, Integer> w1, Map.Entry<String, Integer> w2) -> {
				return w1.getKey().compareTo(w2.getKey());
			});
		});
		
		b2.setOnAction(event -> {
			words.sort((Map.Entry<String,Integer> w1, Map.Entry<String,Integer> w2) -> {
				if (w1.getValue() < w2.getValue()) {
					return 1;
				} else if (w1.getValue() > w2.getValue()) {
					return -1;
				} else {
					return 0;
				}
			});
		});
		
		b3.setOnAction(event -> {
			String search = tf.getText().trim();
			for (Map.Entry<String, Integer> m : words) {
				if (search.equalsIgnoreCase(m.getKey())) {
					listView.scrollTo(m);
					SelectionModel<Map.Entry<String, Integer>> sm = listView.getSelectionModel();
					sm.clearSelection();
					sm.select(m);
				}
			}
		});
		
		hbox.getChildren().addAll(b1, b2, tf, b3);
		root.setBottom(hbox);
		
		primaryStage.setTitle("BookReader");
		primaryStage.setScene(scene);
		primaryStage.show();
	}

	public static void main(String[] args) {
		Application.launch(args);
	}

}
