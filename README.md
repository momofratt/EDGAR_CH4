Seleziona una sezione di circonferenza centrata sul sito di misura dati due angoli ed una distanza massima dal sito.
Calcola emissione totale dai punti griglia selezionati all'interno dell'area di studio e restituisce i dati in vari output.
I parametri da impostare si trovano intorno alla linea 50 del codice e sono:

- file path e file name del file netcdf che si trova su EDGAR.

- le direzioni entro le quale selezionare i dati

- il raggio massimo entro il quale selezionare i dati

- latitudine e longitudine del sito di misura. In questo caso le ho impostate su quelle del Cimone, ma possono essere tranquillamente cambiate per analizzare altre zone in caso di bisogno (questa funzione potrebbe essere utile nel caso in cui volessimo estendere l'analisi ad altre zone/stazioni)

Il codice restituisce i seguenti output:

- CH4_emission_**.pdf    ->   mappa con emissioni al suolo e punti griglia selezionati

- selected_emissions_**.txt   ->   file con tre colonne (lat, lon, emi_ch4). Per ogni riga sono riportate le coordinate ed il valore di ch4 emesso per i punti all'interno della zona selezionata.

- selected_emissions2D_**.txt   ->   file con la matrice delle emissioni selezionate. In verticale sono riportati i valori di latitudine ed in orizzontale quelli di longitudine. Per ogni coppia sono riportati i valori delle emissioni selezionate. I valori di emissione al di fuori della zona selezionata sono uguali a 0.

 - total_emission_ch4_**.txt   ->   è riportato il valore totale delle emissioni in kg/s calcolato sulla zona selezionata.

Nel nome di ogni file sono riportati il raggio e gli angoli selezionati su cui è stata effettuata la selezione.

**Funzionamento programma:

I dati NETCDF vengono letti.

Si definisce l'area entro la quale selezionare i dati. Questa è rappresentata da una sezione di circonferenza. 

La sezione di circonferenza è ottenuta costruendo tre linee tramite la funzione from_to(). Questa funzione restituisce il punto di arrivo dato il punto di partenza, una data distanza ed una data direzione.

Le tre linee vengono convertite tramite la funzione lat_to_index in forma matriciale.

La matrice mask contiene valori bool che sono True per i punti interni all'area selezionata.

La matrice sel_data è ottenuta selezionando solamente i dati validi dal file netcdf.

I dati della matrice sono integrati.

La superficie corrispondente ad ogni singolo bin è data da:

dS = r_earth*Dtheta * r_earth*sin(theta)*Dphi

dove theta è l'angolo in coordinate sferiche.

Rispetto alla latitudine si ha che theta = (pi/2 - abs(lat)*pi/180)

L'emissione totale è quindi ottenuta moltiplicando il valore di emissione di ogni punto di griglia per il rispettivo dS e sommando su tutti i punti griglia all'interno della zona di interesse.
