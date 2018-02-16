create table mutations (id integer primary key auto_increment, name varchar(50) not null);
create table biclusters (id integer primary key auto_increment, mutation_id integer not null references mutations, name varchar(50) not null);
create table genes (id integer primary key auto_increment, name varchar(80) not null);
create table tfs (id integer primary key auto_increment, name varchar(50) not null);
create table bc_mutation_tf_roles (id integer primary key auto_increment, name varchar(50) not null);
create table bc_tf_bc_roles (id integer primary key auto_increment, name varchar(50) not null);

create table bicluster_genes (bicluster_id integer not null references biclusters, gene_id integer not null references genes);

create table bc_mutation_tf (id integer primary key auto_increment, bicluster_id integer not null references biclusters, mutation_id integer not null references mutations, tf_id integer not null references tfs, role integer not null references bc_mutation_tf_roles);

create table bc_tf (id integer primary key auto_increment, bicluster_id integer not null references biclusters, tf_id integer not null references tfs, role integer not null references bc_tf_roles);

insert into bc_mutation_tf_roles (id, name) values (1, 'down-regulates');
insert into bc_mutation_tf_roles (id, name) values (2, 'up-regulates');

insert into bc_tf_roles (id, name) values (1, 'activates');
insert into bc_tf_roles (id, name) values (2, 'represses');
