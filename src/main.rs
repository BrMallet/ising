use rand::Rng;
use std::borrow::Borrow;

use plotters::prelude::*;
//use std::error::Error;

use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

#[derive(Debug, PartialEq)]

// Criaçao do spin

enum SpinValue {
    Positive,
    Negative,
}

//--------------------

// Flip de Spin

fn negate_spin_value(v: &SpinValue) -> SpinValue {
    if *v == SpinValue::Positive {
        SpinValue::Negative
    } else {
        SpinValue::Positive
    }
}

//--------------------------------

// Criação da rede 2d aleatoria

fn random_2d_network(rows: usize, columns: usize) -> Vec<Vec<SpinValue>> {
    let mut rng = rand::thread_rng();

    let mut spins: Vec<Vec<SpinValue>> = Vec::new();
    // Vec::with_capacity(rows);
    for i in 0..rows {
        spins.push(Vec::new());
        for _ in 0..columns {
            let n1: bool = rng.gen();
            spins[i].push( if n1 { SpinValue::Positive } else { SpinValue::Negative } );
        }
    }
    spins
}

//------------------------------

// Inversão de Spin

fn should_change(delta: f64, t: f64) -> bool {

    if t == 0.0 {false}
    else {
    let mut rng = rand::thread_rng();
    let change = rng.gen_range(0 as f64, 1 as f64);
    if change < (-delta / t).exp() { true } else { false }}
}

//--------------------------



// Termos vizinhos

fn next_index(current_index: usize, max_index: usize) -> usize {
    if current_index == max_index - 1 {
        0
    } else {
        current_index + 1
    }
}

fn prev_index(current_index: usize, max_index: usize) -> usize {
    if current_index == 0 {
        max_index - 1
    } else {
        current_index - 1
    }
}

//------------------

//Spin Total novo 

fn s_num(s: &SpinValue) -> f64 {


    if *s == SpinValue::Positive {
        let h_1: f64 = -1.0;  
        h_1    
    } else {
        let h_1: f64 = 1.0;        
        h_1
    }
    
}

//----------------------------

//Interaçao visinha nova

fn visinho_num (s: &SpinValue,s_visinho: &SpinValue) -> f64 {
    if *s == *s_visinho {
        let h_1: f64 = -1.0;  
        h_1
    } else {
        let h_1: f64 = 1.0;  
        h_1 
    }
}

//-------------------------

//Main

fn main() {   

    // Montagem da rede

    let net_size_r = 40; let net_size_c = 40;   
    let root = BitMapBackend::gif("grafico/3.gif", (800, 600), 1_000).unwrap().into_drawing_area();
    

    //--------------------------------

    //Parametros fisicos
    
    let j: f64 = 1.0; let h: f64 = 0.0;       
    
    //Vetores de data
         
     let mut data_h_medio: Vec<Vec<f64>> = Vec::new();
     let mut data_m: Vec<Vec<f64>> =  Vec::new();
     let mut data_c_v: Vec<Vec<f64>> = Vec::new();
     let mut data_x_m: Vec<Vec<f64>> = Vec::new();

    //-------------------------------
    
    //Representação grafica dos spins

    for c in 0 .. 50 {  

        root.fill(&WHITE).unwrap(); 
        let root = root.margin(10, 10, 10, 10);
        

        let mut chart = ChartBuilder::on(&root)
            .caption(
                format!("Temperatura = {} ", c*10),
                ("sans-serif", 10),
            )
            .build_cartesian_2d(0.0..40.0, 0.0..40.0).unwrap();
       
        


        let mut random_network = random_2d_network(net_size_r, net_size_c);

            for i in 0 .. 40 {
                for j in 0 .. 40{

                    let i_ind = i;
                    let j_ind = j;

                    let s_ij = random_network[i_ind][j_ind].borrow();
                    
                    let k = i as f64;
                    let l = j as f64;
    
                    chart.draw_series(PointSeries::of_element(
                        vec![(k , l)],
                        2,
                        if *s_ij == SpinValue::Positive {                
                        &RED
                        }
                        else {
                            &GREEN
                        },
                        &|c, s, st| {
                            return EmptyElement::at(c)    
                            + Circle::new((0,0),s,st.filled()) 
                            
                        },
                    )).unwrap();
                        
                }
            }
            root.present().unwrap();
    //-----------------------------------------------

    //Montagems dos vetores de data

        let mut t = c as f64 ;
        t = t/10.0 ;
 

        

        let mut data_t1: Vec<f64> = Vec::new();
        let mut data_t2: Vec<f64> = Vec::new();
        let mut data_t3: Vec<f64> = Vec::new();
        let mut data_t4: Vec<f64> = Vec::new();   

        data_t1.push(t) ;
        data_t2.push(t) ;
        data_t3.push(t) ;
        data_t4.push(t) ;        

        data_m.push(data_t1);
        data_h_medio.push(data_t2);  
        data_c_v.push(data_t3);        
        data_x_m.push(data_t4);

    // -----------------------------
    
    //Loop -----------------------

        for z in 0 .. 100000000   {        
            
     

        //Spin aleatorio e primeiros vizinhos

        let mut rng = rand::thread_rng();


        let i_ind: usize = rng.gen_range(0 as usize, net_size_r);
        let j_ind: usize = rng.gen_range(0 as usize, net_size_c); 
        let s_ij = random_network[i_ind][j_ind].borrow();

        let i_ind_plus1 = next_index(i_ind, net_size_r);
        let i_ind_minus1 = prev_index(i_ind, net_size_r);
        let j_ind_plus1 = next_index(j_ind, net_size_c);
        let j_ind_minus1 = prev_index(j_ind, net_size_c);

        let s_iplus1_j = random_network[i_ind_plus1][j_ind].borrow();
        let s_iminus1_j = random_network[i_ind_minus1][j_ind].borrow();
        let s_i_jplus1 = random_network[i_ind][j_ind_plus1].borrow();
        let s_ijminus1 = random_network[i_ind][j_ind_minus1].borrow();
        
        //-------------------------

        //Calculo de Delta

        let h_1 = s_num(s_ij)*2.0*h;
        let h_2 = visinho_num(s_ij, s_iplus1_j)*2.0*j;
        let h_3 = visinho_num(s_ij, s_iminus1_j)*2.0*j;
        let h_4 = visinho_num(s_ij, s_i_jplus1)*2.0*j;
        let h_5 = visinho_num(s_ij, s_ijminus1)*2.0*j;

        let delta: f64 = -(h_1 + h_2 + h_3 + h_4 + h_5) ;

        //--------------------------

        //Inversão do Spin

        if delta <= 0.0 {
            random_network[i_ind][j_ind] = negate_spin_value(random_network[i_ind][j_ind].borrow())
        } else {
            if should_change(delta, t) {
                random_network[i_ind][j_ind] = negate_spin_value(random_network[i_ind][j_ind].borrow())
            }
        }

        //Representação grafica em 3 intervalos

        if z == 25000000 {
            for i in 0 .. 40 {
                for j in 0 .. 40{

                    let i_ind = i;
                    let j_ind = j;

                    let s_ij = random_network[i_ind][j_ind].borrow();
                    
                    let k = i as f64;
                    let l = j as f64;
    
                    chart.draw_series(PointSeries::of_element(
                        vec![(k , l)],
                        2,
                        if *s_ij == SpinValue::Positive {                
                        &RED
                        }
                        else {
                            &GREEN
                        },
                        &|c, s, st| {
                            return EmptyElement::at(c)    
                            + Circle::new((0,0),s,st.filled()) 
                            
                        },
                    )).unwrap();
                        
                }
            }
            root.present().unwrap();

        }

        if z == 50000000 {
            for i in 0 .. 40 {
                for j in 0 .. 40{

                    let i_ind = i;
                    let j_ind = j;

                    let s_ij = random_network[i_ind][j_ind].borrow();
                    
                    let k = i as f64;
                    let l = j as f64;
    
                    chart.draw_series(PointSeries::of_element(
                        vec![(k , l)],
                        2,
                        if *s_ij == SpinValue::Positive {                
                        &RED
                        }
                        else {
                            &GREEN
                        },
                        &|c, s, st| {
                            return EmptyElement::at(c)   
                            + Circle::new((0,0),s,st.filled()) 
                            
                        },
                    )).unwrap();
                        
                }
            }
            root.present().unwrap();
        }
        if z == 75000000{
            for i in 0 .. 40 {
                for j in 0 .. 40{

                    let i_ind = i;
                    let j_ind = j;

                    let s_ij = random_network[i_ind][j_ind].borrow();
                    
                    let k = i as f64;
                    let l = j as f64;
    
                    chart.draw_series(PointSeries::of_element(
                        vec![(k , l)],
                        2,
                        if *s_ij == SpinValue::Positive {                
                        &RED
                        }
                        else {
                            &GREEN
                        },
                        &|c, s, st| {
                            return EmptyElement::at(c)    
                            + Circle::new((0,0),s,st.filled()) 
                            
                        },
                    )).unwrap();
                        
                }
            }
            root.present().unwrap();
        } 
        //--------------------------------------------------  
           
        } 

    //------------------------

    //Valores Medios das Grandezas Termodinamicas

        let mut s_total: f64 = 0.0; let mut alpha: f64 = 0.0; let mut h_sqr: f64 = 0.0;

        for l in 0..net_size_r {
            for k in 0..net_size_c {
                let i_ind = l;
                let j_ind = k;

                let s_ij = random_network[i_ind][j_ind].borrow();

                //Calculo númerico
                
                let s = s_num(s_ij);
                s_total = s_total + s ;

                let i_ind_plus1 = next_index(i_ind, net_size_r);
                let i_ind_minus1 = prev_index(i_ind, net_size_r);
                let j_ind_plus1 = next_index(j_ind, net_size_c);
                let j_ind_minus1 = prev_index(j_ind, net_size_c);

                let s_iplus1_j = random_network[i_ind_plus1][j_ind].borrow();
                let s_iminus1_j = random_network[i_ind_minus1][j_ind].borrow();
                let s_i_jplus1 = random_network[i_ind][j_ind_plus1].borrow();
                let s_ijminus1 = random_network[i_ind][j_ind_minus1].borrow();
                
                let h_2 = visinho_num(s_ij, s_iplus1_j);
                let h_3 = visinho_num(s_ij, s_iminus1_j);
                let h_4 = visinho_num(s_ij, s_i_jplus1);
                let h_5 = visinho_num(s_ij, s_ijminus1);

                let vizinho = h_2 + h_3 + h_4 + h_5 ;

                alpha = (alpha + vizinho)/2.0 ;

                let sqr : f64 = (-j * vizinho - h * s).powf(2.0);

                h_sqr  = sqr + h_sqr

            }
        }
    
        let n_c: f64 = net_size_c as f64;
        let n_r: f64 = net_size_r as f64;

        //Magnetização média

        let m: f64 = s_total/(n_c*n_r);        

        data_m[c].push(m) ;

        //Energia livre média

        let h_medio: f64 = (-j * alpha - h * s_total)/(n_c*n_r);

        data_h_medio[c].push(h_medio) ;

        //Capacidade térmica

        let c_v: f64 = (h_sqr/(n_c*n_r) - (h_medio).powf(2.0))/(t).powf(2.0);
        
        data_c_v[c].push(c_v) ;

        //Suscetibilidade magnética

        let x_m: f64 = (n_c*n_r - (m).powf(2.0))/(t).powf(2.0) ;

        data_x_m[c].push(x_m) ;

        

       /* println!("({:?} , {:?} ),", t, m);   
        
        println!("({:?} , {:?} )", t, h_medio );  

        println!("({:?} , {:?} ),", t, c_v);        
        
        println!("({:?} , {:?} ),", t, x_m ); */ 
        
        //Criação dos arquivos dos dados

        let m_str = format!("{:?}", data_m);

        let h_medio_str = format!("{:?}", data_h_medio);

        let c_v_str = format!("{:?}", data_c_v);

        let x_m_str = format!("{:?}", data_x_m);

        let path = Path::new("data_m.txt");
        let display = path.display();
   
       
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => file,
        };
   
       
        match file.write_all(m_str.as_bytes()) {
           Err(why) => panic!("couldn't write to {}: {}", display, why),
           Ok(_) => println!("successfully wrote to {}", display),
       } 

       let path = Path::new("data_h_medio.txt");
       let display = path.display();
   
       
       let mut file = match File::create(&path) {
           Err(why) => panic!("couldn't create {}: {}", display, why),
           Ok(file) => file,
       };
   
       
       match file.write_all(h_medio_str.as_bytes()) {
           Err(why) => panic!("couldn't write to {}: {}", display, why),
           Ok(_) => println!("successfully wrote to {}", display),
       } 

       let path = Path::new("data_c_v.txt");
       let display = path.display();
   
       
       let mut file = match File::create(&path) {
           Err(why) => panic!("couldn't create {}: {}", display, why),
           Ok(file) => file,
       };
   
       
       match file.write_all(c_v_str.as_bytes()) {
           Err(why) => panic!("couldn't write to {}: {}", display, why),
           Ok(_) => println!("successfully wrote to {}", display),
       } 


       let path = Path::new("data_x_m.txt");
       let display = path.display();
   
       
       let mut file = match File::create(&path) {
           Err(why) => panic!("couldn't create {}: {}", display, why),
           Ok(file) => file,
       };
   
       
       match file.write_all(x_m_str.as_bytes()) {
           Err(why) => panic!("couldn't write to {}: {}", display, why),
           Ok(_) => println!("successfully wrote to {}", display),
       }

       //-----------------------------------------------

       //Ultima representação gráfica

       for i in 0 .. 40 {
        for j in 0 .. 40{

            let i_ind = i;
            let j_ind = j;

            let s_ij = random_network[i_ind][j_ind].borrow();
            
            let k = i as f64;
            let l = j as f64;

            chart.draw_series(PointSeries::of_element(
                vec![(k , l)],
                2,
                if *s_ij == SpinValue::Positive {                
                &RED
                }
                else {
                    &GREEN
                },
                &|c, s, st| {
                    return EmptyElement::at(c)    
                    + Circle::new((0,0),s,st.filled()) 
                    
                },
            )).unwrap();
                
        }
    }
    root.present().unwrap(); 
    //---------------------------------------------------

    


        
        
    } 

   



    
    


      //-------------------------------------------------      

}