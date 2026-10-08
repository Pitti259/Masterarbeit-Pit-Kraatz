# Skalierbarkeit der Satellitenfertigung

Vereinfachtes Ablaufmodell zur Darstellung von Zusammenhängen der verschiedenen Produktionsmerkmalen (Kapitel 6.3 der Masterarbeit). 

Das Repository macht die Ergebnisse nachvollziehbar. Es enthält den
vollständigen Quellcode, die Ergebnisse eines Laufs und eine aufbereitete Ergebnisübersicht.

## Inhalt

- `Simulation - Satellitenfertigung.py` Modell, Verifikation, vier Versuchsreihen, Sensitivitätsanalyse
- `Simulationsprotokoll.txt` Ausgabe eines vollständigen Laufs
- `Übersicht der gesamten Ergebnisse.xlsx` alle Ergebnisse aufbereitet, mit Lesehilfe
- `Ergebnisse - *.csv` und `Sensitivitätsanalyse.csv` die Rohwerte je Versuchsreihe

## Nutzung

Zum Nachlesen der Ergebnisse eignet sich die Excel-Übersicht, das erste Blatt erklärt Aufbau und Begriffe. 
Einen schnellen Überblick gibt das Simulationsprotokoll. 
Die CSV-Dateien enthalten dieselben Zahlen unformatiert und mit allen Nachkommastellen.

Zum Ausführen des Modells:

    python "Simulation - Satellitenfertigung.py"

Python 3.9 oder neuer, keine zusätzlichen Pakete, Laufzeit rund eineinhalb
Minuten. Modellannahmen und Parameterwerte stehen in Anhang F der Arbeit.
