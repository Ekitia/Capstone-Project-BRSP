# Study Case_OmicsLite_Capstone Project

## Proses Pengerjaan Analisis Data DEG dari GEO NCBI
<img alt="image" src="https://github.com/Ekitia/Capstone-Project-BRSP/blob/main/Sumber%20Gambar/unnamed.jpg" />

## Topic Capstone Project: Kondisi Gen Perokok
### Tools: GEO, GEO2R, R (R/Rstudio), GO (g:profiler), KEGG 
Proses tersebut dilakukan dalam mengolah data DEG pada GEO NCBI (GSE40885) dimana hasil akhir diinpretasikan berdasarkan pelatihan yang dipelajari pada Omicslite

### Grup yang Digunakan dalam Penelitian: Saline (Sal) dan LPS

## Langkah-Langkah Analisis GEO / GEO2R
1. Melakukan pengambilan data DEG melalui website GEO NCBI dengan kasus utama berupa gen perokok akut [Link GEO2R](<https://www.ncbi.nlm.nih.gov/geo/geo2r/>)
3. Dilakukan pemeriksaan dataset (GSE40885), lalu dilakukan analisis data dengan GEO2R di NCBI
4. Di situs web analisis GEO2R, ​​dapat dilakukan filter baris dan pemberian nama, yaitu pada lokasi pada cairan tubuh (Saline) dan Lipopolisakarida (LPS)
5. Dilakukan pemeriksaan opsi analisis di taskbar (di bawah)
<img alt="image" src=https://github.com/Ekitia/Capstone-Project-BRSP/blob/main/Sumber%20Gambar/543415513-8a183f21-f064-4201-a463-98cd662eb183.png />
6. Pilih "Opsi", Biasanya secara default opsi sudah dipilih. Tetapi, terapkan "Disediakan oleh Pengirim" di kategori platform (Karena NCBI mungkin menghasilkan beberapa perubahan pada dataset) dan pilih "ya" untuk menerapkan limma.
|
Hasil Akhir:

<img alt="image" src="https://github.com/Ekitia/Capstone-Project-BRSP/blob/main/Sumber%20Gambar/Screenshot%20from%202026-03-08%2021-49-19.png" />
|
## Langkah-Langkah Analisis GEO dengan RStudio
1. Melakukan pengambilan data GEO menggunakan package GEOquery sesuai nomor GEO yang telah ditentukan
2. Membuat grup data sampel berdasarkan lokasi di Saline dan LPS
3. Membuat sistem matriks data Sampel
4. Melakukan Analisis Ekspresi Differensial (LIMMA)
5. Memberi ulang nama gen dari sampel yang telah dianalisis LIMMA
6. Visualisasi Data DEG (barplot, volcano plot, UMAP, dan lain-lain), serta menyimpan hasil dataset dalam format .csv
|
Hasil Akhir:
<img alt="image" src="https://github.com/Ekitia/Capstone-Project-BRSP/blob/main/Sumber%20Gambar/Screenshot%20from%202026-03-08%2021-52-17.png" />
|
## Langkah Tambahan untuk Kedua Proses Data Tersebut
1. Persiapan dan Penyaringan Data (Data Preparation)
a) Penyaringan DEG: Dari dataset awal, saring 20 gen yang mengalami peningkatan ekspresi (upregulated) dan 20 gen yang mengalami penurunan ekspresi (downregulated). Anda bisa menggunakan bantuan alat seperti ChatGPT (dengan prompt yang sesuai), SQL, atau Microsoft Excel untuk memudahkan konversi (.tsv ke .csv) dan pemfilteran.
b) Dibersihkan karakter tersebut (misalnya menggunakan ChatGPT) dan susun daftar gen secara vertikal agar siap dianalisis.

3. Analisis Gene Ontology (GO)
a) Input Data: Kunjungi situs web alat analisis GO seperti g:Profiler atau DAVID, lalu masukkan daftar gen vertikal yang sudah disiapkan.
b) Eksekusi dan Eksplorasi: Jalankan analisis dan periksa hasilnya pada kategori GO-BP (Biological Process), GO-MP (Molecular Function), dan GO-CC (Cellular Component).
c) Identifikasi Proses: Cari proses biologis utama yang relevan dengan sampel Anda, seperti "Respons imun bawaan", "Respons inflamasi", "Respons terhadap virus", atau "Proses apoptosis", dan eksplorasi hasil lainnya untuk dasar kesimpulan.

4. Pemetaan Jalur KEGG (KEGG Pathway Analysis)
a) Input Terpisah: Buka situs KEGG Mapper dan atur pencarian untuk organisme Homo sapiens. Masukkan data gen yang mengalami penurunan ekspresi terlebih dahulu untuk dianalisis, lalu ulangi proses yang sama secara terpisah untuk gen yang mengalami peningkatan ekspresi.
b) Evaluasi Jalur: Analisis output jalur yang muncul (misalnya: "Jalur pensinyalan reseptor Toll-like", "Interaksi sitokin-reseptor sitokin", atau "Apoptosis") beserta jumlah gen yang ditemukan di dalamnya.
Integrasi: Gabungkan temuan jalur KEGG ini dengan hasil analisis GO sebelumnya.

5. Interpretasi dan Kesimpulan Akhir
a) Pengumpulan Visualisasi: Kumpulkan semua data visualisasi yang Anda miliki, yang mencakup: Data visualisasi awal dari GEO2R (menggunakan R/Rstudio), daftar hasil pengayaan dari GO (BP, MP, CC), serta data visualisasi jalur gen dari KEGG.
b) Penarikan Kesimpulan: Gunakan seluruh kumpulan data dan visualisasi tersebut secara holistik untuk menyusun kesimpulan akhir dari penelitian atau analisis Anda.
