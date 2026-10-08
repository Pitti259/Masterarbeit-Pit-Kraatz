# Ereignisdiskretes Ablaufmodell zur Skalierbarkeit der Satellitenfertigung
# Demonstrator zu Kapitel 6.3
#
# Saemtliche Modellannahmen und Parameterwerte sind in Anhang F dokumentiert
# und werden hier nicht wiederholt. Alle Zeiten sind normierte Zeiteinheiten (ZE).
#
# Ausfuehrung: python simulation_satellitenfertigung.py
# Es werden keine zusaetzlichen Pakete benoetigt, Laufzeit rund eineinhalb Minuten.
# Ergebnisse werden als protokoll.txt und als CSV-Dateien im selben Ordner abgelegt.

import csv
import heapq
import math
import os
import random
import statistics
from itertools import count


# -------------------------------------------------------------
# EREIGNISDISKRETER KERN
# -------------------------------------------------------------

# Zeitpunkt, an dem wartende Prozesse fortgesetzt werden
class Ereignis:
    __slots__ = ("umgebung", "rueckrufe", "eingetreten")

    def __init__(self, umgebung):
        self.umgebung = umgebung
        self.rueckrufe = []
        self.eingetreten = False


# Terminkalender und Simulationsuhr
class Umgebung:

    def __init__(self):
        self.jetzt = 0.0
        self._kalender = []
        self._zaehler = count()

    def eintragen(self, ereignis, verzoegerung=0.0):
        heapq.heappush(self._kalender,
                       (self.jetzt + verzoegerung, next(self._zaehler), ereignis))

    def warten(self, dauer):
        ereignis = Ereignis(self)
        self.eintragen(ereignis, dauer)
        return ereignis

    def prozess(self, ablauf):
        return Prozess(self, ablauf)

    # Ereignisse der Reihe nach abarbeiten, bis der Kalender leer ist
    def laufen(self):
        while self._kalender:
            zeitpunkt, _, ereignis = heapq.heappop(self._kalender)
            self.jetzt = zeitpunkt
            ereignis.eingetreten = True
            rueckrufe, ereignis.rueckrufe = ereignis.rueckrufe, []
            for rueckruf in rueckrufe:
                rueckruf(ereignis)


# Ablauf, der an Ereignissen unterbrochen und fortgesetzt wird
class Prozess(Ereignis):

    def __init__(self, umgebung, ablauf):
        super().__init__(umgebung)
        self.ablauf = ablauf
        start = Ereignis(umgebung)
        start.rueckrufe.append(self._fortsetzen)
        umgebung.eintragen(start)

    def _fortsetzen(self, _):
        try:
            naechstes = next(self.ablauf)
        except StopIteration:
            return
        naechstes.rueckrufe.append(self._fortsetzen)


# Belegbare Einheit mit Warteschlange und Auslastungserfassung
class Ressource:

    def __init__(self, umgebung, kapazitaet=1):
        self.umgebung = umgebung
        self.kapazitaet = kapazitaet
        self.belegt = 0
        self.warteschlange = []
        self._flaeche = 0.0
        self._letzte = 0.0

    def _fortschreiben(self):
        self._flaeche += self.belegt * (self.umgebung.jetzt - self._letzte)
        self._letzte = self.umgebung.jetzt

    def belegen(self):
        ereignis = Ereignis(self.umgebung)
        if self.belegt < self.kapazitaet:
            self._fortschreiben()
            self.belegt += 1
            self.umgebung.eintragen(ereignis)
        else:
            self.warteschlange.append(ereignis)
        return ereignis

    def freigeben(self, _):
        self._fortschreiben()
        self.belegt -= 1
        if self.warteschlange:
            naechstes = self.warteschlange.pop(0)
            self._fortschreiben()
            self.belegt += 1
            self.umgebung.eintragen(naechstes)

    def auslastung(self, gesamtdauer):
        self._fortschreiben()
        if gesamtdauer <= 0:
            return 0.0
        return self._flaeche / (self.kapazitaet * gesamtdauer)


# -------------------------------------------------------------
# MODELLPARAMETER
# -------------------------------------------------------------
# Herkunft, Variation und Begruendung der Werte: siehe Anhang F

P = dict(
    # Programm und Ablauf
    N=180,
    stationsplaetze=30,
    W_0=100.0,

    # Produktionsmerkmale
    f_tiefe=0.65,
    k_auto=0.80,
    r_ruest=0.20,
    v_wechsel=0.30,

    # Streuung
    cv=0.30,

    # Lernkurve
    lernrate=0.85,
    t_min_anteil=0.75,

    # Zulieferung
    module=6,
    p_fehlteil=0.30,
    t_zuliefer=3.0,

    # Gemeinsame Pruefressource
    pruefstationen=3,
    t_pruef=2.5,

    # Durchfuehrung
    wiederholungen=30,
    startwert=20260906,
)

KONFIGURATIONEN = [(1, 30), (2, 15), (3, 10), (5, 6), (6, 5), (10, 3), (15, 2)]


# -------------------------------------------------------------
# HILFSFUNKTIONEN
# -------------------------------------------------------------

# Lernrate in den Exponenten b umrechnen
def lernexponent(lernrate):
    return -math.log(lernrate) / math.log(2.0)


# Lognormalverteilte Zeit mit gegebenem Mittelwert und Variationskoeffizienten
def lognormal(zufall, mittelwert, cv):
    if cv <= 0 or mittelwert <= 0:
        return max(mittelwert, 0.0)
    sigma = math.sqrt(math.log(1 + cv * cv))
    return zufall.lognormvariate(math.log(mittelwert) - 0.5 * sigma * sigma, sigma)


# -------------------------------------------------------------
# MODELL EINES PROGRAMMDURCHLAUFS
# -------------------------------------------------------------

# Simuliert ein Programm mit L parallelen Linien zu je S Stationen.
# Ueber abweichung lassen sich einzelne Parameter fuer ein Experiment
# ueberschreiben, zum Beispiel simuliere(5, 6, k_auto=0.6).
def simuliere(L, S, p=P, startwert=0, **abweichung):
    w = dict(p)
    w.update(abweichung)
    zufall = random.Random(startwert)
    umgebung = Umgebung()

    # Schritt 1: Produktionsmerkmale in Zeitgroessen ueberfuehren
    W_int = w["W_0"] * w["f_tiefe"]
    t_nom = W_int / S
    t_0 = w["k_auto"] * t_nom
    t_min = w["t_min_anteil"] * t_0
    t_ruest = w["r_ruest"] * t_nom
    b = lernexponent(w["lernrate"])
    p_fehl = (1.0 - w["f_tiefe"]) * w["p_fehlteil"]

    # Schritt 2: Ressourcen und Zaehler anlegen
    stationen = [[Ressource(umgebung) for _ in range(S)] for _ in range(L)]
    pruefung = Ressource(umgebung, w["pruefstationen"])
    gefertigt = [[0] * S for _ in range(L)]
    offen = list(range(w["N"]))
    wechsel = {}
    fehlteile = {}
    fertig = []
    zaehler_ruest = [0]
    zaehler_fehl = [0]

    # Schritt 3: Zeitbedarf eines Satelliten an einer Station
    def stationszeit(linie, station, satellit):
        gefertigt[linie][station] += 1
        n = gefertigt[linie][station]
        zeit = t_min + (t_0 - t_min) * n ** (-b)          # Lernkurve
        if wechsel[satellit]:                             # Ruesten
            zeit += t_ruest
            if station == 0:
                zaehler_ruest[0] += 1
        zeit = lognormal(zufall, zeit, w["cv"])           # Streuung
        for _ in range(fehlteile[satellit][station]):     # Warten auf Fehlteil
            zeit += zufall.expovariate(1.0 / w["t_zuliefer"])
        return zeit

    # Schritt 4: Durchlauf eines Satelliten durch seine Linie und die Pruefung
    def durchlauf(linie, satellit, belegung):
        for station in range(S):
            yield umgebung.warten(stationszeit(linie, station, satellit))
            if station + 1 < S:
                # ohne Puffer: erst die naechste Station belegen, dann freigeben
                naechste = stationen[linie][station + 1].belegen()
                yield naechste
                stationen[linie][station].freigeben(belegung)
                belegung = naechste
            else:
                stationen[linie][station].freigeben(belegung)
        platz = pruefung.belegen()
        yield platz
        yield umgebung.warten(lognormal(zufall, w["t_pruef"], w["cv"]))
        pruefung.freigeben(platz)
        fertig.append(umgebung.jetzt)

    # Schritt 5: Satelliten an die Linien uebergeben
    def zuteilung(linie):
        erster = True
        while offen:
            belegung = stationen[linie][0].belegen()
            yield belegung
            if not offen:
                stationen[linie][0].freigeben(belegung)
                return
            satellit = offen.pop(0)
            wechsel[satellit] = (not erster) and (zufall.random() < w["v_wechsel"])
            erster = False
            treffer = [0] * S
            for _ in range(w["module"]):
                if zufall.random() < p_fehl:
                    treffer[zufall.randrange(S)] += 1
                    zaehler_fehl[0] += 1
            fehlteile[satellit] = treffer
            umgebung.prozess(durchlauf(linie, satellit, belegung))

    for linie in range(L):
        umgebung.prozess(zuteilung(linie))
    umgebung.laufen()

    # Schritt 6: Ergebnisgroessen bestimmen
    fertig.sort()
    programmdauer = fertig[-1]
    fenster = max(10, w["N"] // 9)
    raten = [fenster / (fertig[i] - fertig[i - fenster])
             for i in range(fenster, len(fertig))]
    endrate = statistics.fmean(raten[-max(3, fenster // 2):])
    ziel = 0.95 * endrate
    hochlauf = programmdauer
    for i, rate in enumerate(raten):
        if rate >= ziel:
            hochlauf = fertig[i + fenster]
            break
    return dict(
        L=L, S=S,
        programmdauer=programmdauer,
        endrate=endrate,
        hochlauf=hochlauf,
        auslastung_pruefung=pruefung.auslastung(programmdauer),
        ruestquote=zaehler_ruest[0] / w["N"],
        fehlteile_je_satellit=zaehler_fehl[0] / w["N"],
        fertiggestellt=len(fertig),
    )


# -------------------------------------------------------------
# WIEDERHOLUNGEN UND KONFIDENZINTERVALLE
# -------------------------------------------------------------

KENNZAHLEN = ("programmdauer", "endrate", "hochlauf",
              "auslastung_pruefung", "ruestquote", "fehlteile_je_satellit")


# Wiederholt eine Konfiguration und liefert Mittelwerte mit 95-%-Intervall.
# Alle Konfigurationen verwenden denselben Satz von Zufallsstartwerten.
def auswerten(L, S, p=P, **abweichung):
    laeufe = [simuliere(L, S, p, startwert=p["startwert"] + 97 * i, **abweichung)
              for i in range(p["wiederholungen"])]
    n = len(laeufe)
    ergebnis = dict(L=L, S=S, wiederholungen=n)
    for kennzahl in KENNZAHLEN:
        werte = [lauf[kennzahl] for lauf in laeufe]
        ergebnis[kennzahl] = statistics.fmean(werte)
        ergebnis[kennzahl + "_ki"] = 1.96 * statistics.stdev(werte) / math.sqrt(n)
    ergebnis["fertiggestellt"] = laeufe[0]["fertiggestellt"]
    return ergebnis


# -------------------------------------------------------------
# AUSGABE
# -------------------------------------------------------------

ORDNER = os.path.dirname(os.path.abspath(__file__))
_protokoll = None


# Ausgabe gleichzeitig in Konsole und Protokolldatei
def notiere(text=""):
    print(text)
    if _protokoll is not None:
        _protokoll.write(text + "\n")


# Ergebniszeilen als CSV ablegen
def speichern(name, zeilen):
    with open(os.path.join(ORDNER, name), "w", newline="") as datei:
        schreiber = csv.DictWriter(datei, fieldnames=list(zeilen[0].keys()))
        schreiber.writeheader()
        schreiber.writerows(zeilen)


# -------------------------------------------------------------
# VERIFIKATION
# -------------------------------------------------------------

# Grenzwerttests, Invarianzpruefungen und analytische Untergrenzen
def verifikation():

    def mittel(L, S, laeufe=30, **abweichung):
        werte = [simuliere(L, S, P, startwert=P["startwert"] + 97 * i,
                           **abweichung)["programmdauer"] for i in range(laeufe)]
        return sum(werte) / len(werte)

    notiere("V  VERIFIKATION UND PLAUSIBILISIERUNG")
    notiere()

    # Einzelne Einfluesse entfernen oder verstaerken und die Richtung pruefen
    notiere("   Grenzwerttests (Konfiguration 5 x 6, 3 Pruefstationen,"
            " je 30 Laeufe, Programmdauer in ZE)")
    notiere(f"     Basislauf                      {mittel(5, 6):7.1f}")
    notiere(f"     ohne Lerneffekt                {mittel(5, 6, lernrate=1.0):7.1f}"
            "   erwartet: laenger")
    notiere(f"     ohne Mindestzeit               {mittel(5, 6, t_min_anteil=0.0):7.1f}"
            "   erwartet: kuerzer")
    notiere(f"     ohne Variantenwechsel          {mittel(5, 6, v_wechsel=0.0):7.1f}"
            "   erwartet: kuerzer")
    notiere(f"     ohne Streuung                  {mittel(5, 6, cv=0.0):7.1f}"
            "   erwartet: kuerzer")
    notiere(f"     volle Fertigungstiefe          {mittel(5, 6, f_tiefe=1.0):7.1f}"
            "   erwartet: laenger")
    notiere(f"     ohne Automatisierungsvorteil   {mittel(5, 6, k_auto=1.0):7.1f}"
            "   erwartet: laenger")
    notiere(f"     eine Pruefstation              {mittel(5, 6, pruefstationen=1):7.1f}"
            "   erwartet: laenger")
    notiere()

    # Groessen, die von der Aufteilung unabhaengig sein muessen
    notiere("   Invarianzen ueber die Konfigurationen (je 30 Laeufe)")
    notiere("     Erwartet werden 0,30 * (180 - L) / 180 Anpassungsvorgaenge und"
            " 0,63 Fehlteile je Satellit.")
    notiere("     Konfig      Anpassungen  erwartet   Fehlteile je Satellit")
    for L, S in KONFIGURATIONEN:
        e = auswerten(L, S)
        soll = P["v_wechsel"] * (P["N"] - L) / P["N"]
        notiere(f"     {L:2d} x {S:2d}       {e['ruestquote']:.3f}       {soll:.3f}"
                f"          {e['fehlteile_je_satellit']:.3f}")
    notiere()

    # Abgleich der Lernkurve mit der Modellrechnung der Quelle
    b = lernexponent(P["lernrate"])
    t648 = 360 + (480 - 360) * 648 ** (-b)
    notiere("   Lernkurve gegen die Quelle (480 / 360 Minuten, 648. Einheit)")
    notiere(f"     b = {b:.4f}   t(648) = {t648:.1f} min   Reduktion "
            f"{(480 - t648) / 480 * 100:.2f} %   (Quelle: 19,4 %)")
    notiere()

    # Rechnerische Untergrenzen, die kein Simulationslauf unterschreiten darf
    W_int = P["W_0"] * P["f_tiefe"]
    bedarf = P["k_auto"] * W_int + P["v_wechsel"] * P["r_ruest"] * W_int
    grenze_fertigung = (P["N"] * (P["t_min_anteil"] * P["k_auto"] * W_int
                                  + P["v_wechsel"] * P["r_ruest"] * W_int)
                        / P["stationsplaetze"])
    notiere("   Analytische Untergrenzen (Referenzfall)")
    notiere(f"     Stationsbedarf je Satellit                   {bedarf:6.1f} ZE")
    notiere(f"     Untergrenze Fertigung bei vollem Lerneffekt  {grenze_fertigung:6.1f} ZE")
    for K in (1, 2, 3):
        notiere(f"     Untergrenze Pruefung bei K = {K}                "
                f"{P['N'] * P['t_pruef'] / K:6.1f} ZE")
    notiere("     Kein simulierter Wert darf diese Grenzen unterschreiten.")
    notiere()


# -------------------------------------------------------------
# EXPERIMENTE E1 BIS E4
# -------------------------------------------------------------

STANDARDISIERUNG = ((0.30, "niedrig"), (0.20, "mittel"), (0.10, "hoch"))
AUTOMATISIERUNG = ((1.00, "niedrig"), (0.80, "mittel"), (0.60, "hoch"))


# E1  Produktionsorganisation x Pruefkapazitaet
def experiment_1():
    notiere("E1  PRODUKTIONSORGANISATION x PRUEFKAPAZITAET")
    zeilen = []
    for K in (1, 2, 3):
        for L, S in KONFIGURATIONEN:
            e = auswerten(L, S, pruefstationen=K)
            e["pruefstationen"] = K
            zeilen.append(e)
            notiere(f"   K = {K}   {L:2d} x {S:2d}   Programmdauer "
                    f"{e['programmdauer']:7.1f} +- {e['programmdauer_ki']:4.1f}"
                    f"   Hochlauf {e['hochlauf']:6.1f}"
                    f"   Auslastung Pruefung {e['auslastung_pruefung'] * 100:5.1f} %")
    speichern("e1_organisation_pruefkapazitaet.csv", zeilen)
    notiere()


# E2  Standardisierung x Automatisierung, je bei knapper und ausreichender Pruefung
def experiment_2():
    notiere("E2  STANDARDISIERUNG x AUTOMATISIERUNG")
    zeilen = []
    for L, S in ((1, 30), (5, 6)):
        for K in (1, 3):
            for r, r_name in STANDARDISIERUNG:
                for k, k_name in AUTOMATISIERUNG:
                    e = auswerten(L, S, r_ruest=r, k_auto=k, pruefstationen=K)
                    e.update(pruefstationen=K, standardisierung=r_name, r_ruest=r,
                             automatisierung=k_name, k_auto=k)
                    zeilen.append(e)
                    notiere(f"   {L:2d} x {S:2d}  K = {K}  Standardisierung "
                            f"{r_name:7s} Automatisierung {k_name:7s}   "
                            f"{e['programmdauer']:7.1f} +- {e['programmdauer_ki']:4.1f}")
    speichern("e2_standardisierung_automatisierung.csv", zeilen)
    notiere()


# E3  Fertigungstiefe, anschliessend Kipppunkt der Lieferverzoegerung
def experiment_3():
    notiere("E3  FERTIGUNGSTIEFE")
    zeilen = []
    for K in (1, 3):
        for f in (0.50, 0.65, 0.80):
            for L, S in KONFIGURATIONEN:
                e = auswerten(L, S, f_tiefe=f, pruefstationen=K)
                e.update(pruefstationen=K, f_tiefe=f)
                zeilen.append(e)
                notiere(f"   K = {K}  f = {f:.2f}   {L:2d} x {S:2d}   "
                        f"{e['programmdauer']:7.1f} +- {e['programmdauer_ki']:4.1f}"
                        f"   Fehlteile je Satellit {e['fehlteile_je_satellit']:.2f}")
    speichern("e3_fertigungstiefe.csv", zeilen)

    notiere()
    notiere("   Kipppunkt: ab welcher Lieferverzoegerung dreht die Rangfolge?"
            " (5 x 6, K = 3)")
    notiere("     Wartezeit    f = 0,50   f = 0,65   f = 0,80")
    zeilen = []
    for t_z in (3, 6, 10, 15, 20, 30):
        werte = {}
        for f in (0.50, 0.65, 0.80):
            e = auswerten(5, 6, f_tiefe=f, t_zuliefer=t_z, pruefstationen=3)
            werte[f] = e["programmdauer"]
            zeilen.append(dict(t_zuliefer=t_z, f_tiefe=f, **e))
        notiere(f"     {t_z:2d} ZE       {werte[0.50]:7.1f}    {werte[0.65]:7.1f}"
                f"    {werte[0.80]:7.1f}")
    speichern("e3b_kipppunkt.csv", zeilen)
    notiere()


# E4  Produktportfolio x Standardisierung
def experiment_4():
    notiere("E4  PRODUKTPORTFOLIO x STANDARDISIERUNG")
    zeilen = []
    for L, S in ((1, 30), (5, 6)):
        for K in (1, 3):
            for v in (0.10, 0.30, 0.50):
                for r, r_name in STANDARDISIERUNG:
                    e = auswerten(L, S, v_wechsel=v, r_ruest=r, pruefstationen=K)
                    e.update(pruefstationen=K, v_wechsel=v, standardisierung=r_name)
                    zeilen.append(e)
                    notiere(f"   {L:2d} x {S:2d}  K = {K}  Variantenwechsel {v:.2f}"
                            f"  Standardisierung {r_name:7s}   "
                            f"{e['programmdauer']:7.1f} +- {e['programmdauer_ki']:4.1f}")
    speichern("e4_portfolio_standardisierung.csv", zeilen)
    notiere()


# -------------------------------------------------------------
# SENSITIVITAETSANALYSE
# -------------------------------------------------------------
# Geprueft wird die Richtung der Aussage, nicht ihre Groesse. Gleichheit gilt
# als bestaetigt, wenn die Differenz kleiner ist als die Summe der beiden
# halben Konfidenzintervalle.

EINSTELLUNGEN = [
    ("Referenzfall", {}),
    ("Mindestzeit 0,65", {"t_min_anteil": 0.65}),
    ("Mindestzeit 0,85", {"t_min_anteil": 0.85}),
    ("Streuung CV 0,15", {"cv": 0.15}),
    ("Streuung CV 0,45", {"cv": 0.45}),
    ("Pruefdauer 1,5 ZE", {"t_pruef": 1.5}),
    ("Pruefdauer 3,5 ZE", {"t_pruef": 3.5}),
]


# Prueft, ob die Befunde von einzelnen gesetzten Annahmen abhaengen
def sensitivitaet():
    notiere("S  SENSITIVITAETSANALYSE")
    notiere()
    notiere("   B1  Die Aufteilung auf mehrere Linien verkuerzt die Programmdauer.")
    notiere("   B2  Der Gewinn der Aufteilung faellt groesser aus, wenn die")
    notiere("       Pruefkapazitaet ausreicht.")
    notiere("   B3  Bei einer einzelnen langen Linie bringt zusaetzliche")
    notiere("       Pruefkapazitaet keinen Vorteil.")
    notiere("   B4  Standardisierung wirkt vor allem bei hohem Variantenanteil.")
    notiere()

    zeilen = []
    for name, abweichung in EINSTELLUNGEN:

        # Vergleichswerte fuer die vier Aussagen erzeugen
        werte = {}
        for K in (1, 3):
            for L, S in ((1, 30), (5, 6), (15, 2)):
                werte[(K, L)] = auswerten(L, S, pruefstationen=K, **abweichung)
        standard = {}
        for v in (0.10, 0.50):
            for r in (0.30, 0.10):
                standard[(v, r)] = auswerten(5, 6, pruefstationen=3, v_wechsel=v,
                                             r_ruest=r, **abweichung)

        def d(K, L):
            return werte[(K, L)]["programmdauer"]

        def ki(K, L):
            return werte[(K, L)]["programmdauer_ki"]

        b1 = d(1, 1) > d(1, 15) and d(3, 1) > d(3, 15)
        gewinn_knapp = d(1, 5) - d(1, 15)
        gewinn_reichlich = d(3, 5) - d(3, 15)
        b2 = gewinn_reichlich > gewinn_knapp
        b3 = abs(d(1, 1) - d(3, 1)) < ki(1, 1) + ki(3, 1)
        wirkung_wenig = (standard[(0.10, 0.30)]["programmdauer"]
                         - standard[(0.10, 0.10)]["programmdauer"])
        wirkung_viel = (standard[(0.50, 0.30)]["programmdauer"]
                        - standard[(0.50, 0.10)]["programmdauer"])
        b4 = wirkung_viel > wirkung_wenig

        # Voraussetzung von B3: die Pruefressource darf nicht ohnehin binden
        pruefbedarf = P["N"] * abweichung.get("t_pruef", P["t_pruef"])
        pruefung_bindet_immer = pruefbedarf > d(3, 1)

        ja = lambda x: "haelt" if x else "HAELT NICHT"
        notiere(f"   --- {name}")
        notiere(f"       1 x 30: K=1 {d(1, 1):6.1f}  K=3 {d(3, 1):6.1f}    "
                f"15 x 2: K=1 {d(1, 15):6.1f}  K=3 {d(3, 15):6.1f}")
        notiere(f"       B1  {ja(b1)}")
        notiere(f"       B2  {ja(b2)}   ({gewinn_knapp:.1f} gegen "
                f"{gewinn_reichlich:.1f} ZE, Faktor "
                f"{gewinn_reichlich / gewinn_knapp:.1f})")
        notiere(f"       B3  {ja(b3)}   (Differenz {abs(d(1, 1) - d(3, 1)):.1f},"
                f" KI-Summe {ki(1, 1) + ki(3, 1):.1f})")
        if pruefung_bindet_immer:
            notiere(f"           Hinweis: Der Pruefbedarf von {pruefbedarf:.0f} ZE"
                    " uebersteigt die Dauer der langen Linie; die Pruefressource")
            notiere("           bindet hier in jeder Konfiguration.")
        notiere(f"       B4  {ja(b4)}   ({wirkung_wenig:.1f} gegen "
                f"{wirkung_viel:.1f} ZE, Faktor {wirkung_viel / wirkung_wenig:.1f})")
        notiere()

        zeilen.append(dict(
            einstellung=name,
            d_1x30_K1=d(1, 1), d_1x30_K1_ki=ki(1, 1),
            d_1x30_K3=d(3, 1), d_1x30_K3_ki=ki(3, 1),
            d_5x6_K1=d(1, 5), d_5x6_K3=d(3, 5),
            d_15x2_K1=d(1, 15), d_15x2_K3=d(3, 15),
            gewinn_K1=gewinn_knapp, gewinn_K3=gewinn_reichlich,
            wirkung_std_v10=wirkung_wenig, wirkung_std_v50=wirkung_viel,
            pruefung_bindet_immer=pruefung_bindet_immer,
            B1=b1, B2=b2, B3=b3, B4=b4))

    speichern("sensitivitaet.csv", zeilen)

    notiere("   Zusammenfassung")
    for schluessel, text in (
            ("B1", "Aufteilung verkuerzt die Programmdauer"),
            ("B2", "Gewinn groesser bei ausreichender Pruefkapazitaet"),
            ("B3", "bei langer Linie ist Pruefkapazitaet wirkungslos"),
            ("B4", "Standardisierung wirkt vor allem bei vielen Varianten")):
        treffer = sum(1 for z in zeilen if z[schluessel])
        notiere(f"     {schluessel}  {text:50s} {treffer} von {len(zeilen)}")
    ausnahmen = [z["einstellung"] for z in zeilen if not z["B3"]]
    if ausnahmen:
        notiere()
        notiere("     B3 gilt nicht, wenn die Pruefressource bereits die langsamste")
        notiere("     Konfiguration bindet: " + ", ".join(ausnahmen))
    notiere()


# -------------------------------------------------------------
# PROGRAMMSTART
# -------------------------------------------------------------

if __name__ == "__main__":
    _protokoll = open(os.path.join(ORDNER, "protokoll.txt"), "w", buffering=1)
    notiere("Simulationsmodell zur Skalierbarkeit der Satellitenfertigung")
    notiere("Alle Zeiten in normierten Zeiteinheiten (ZE).")
    notiere(f"Je Konfiguration {P['wiederholungen']} Wiederholungen;"
            " angegeben sind Mittelwert und halbe Breite des 95-%-Konfidenzintervalls.")
    notiere("=" * 74)
    notiere()
    verifikation()
    experiment_1()
    experiment_2()
    experiment_3()
    experiment_4()
    sensitivitaet()
    notiere("FERTIG")
    _protokoll.close()